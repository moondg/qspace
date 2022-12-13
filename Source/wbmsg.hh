#ifndef __WB_MSG_HH__
#define __WB_MSG_HH__

/* https://stackoverflow.com/questions/24785230/difference-between-msgget-and-mq-open
 * Wb,Oct06,22

  The msgXXX family is inherited from SysV. The mq_XXX family was
  created by POSIX.  Both families are however part of POSIX.

  The two families have somewhat different feature sets the most
  important being that the SysV variant allows for the specification
  of different message types. This provides for flexibility in the
  way messages are processed and in assigning priorities. The POSIX
  family in contrast always returns the oldest message with the
  highest priority.

  As a rule the POSIX family is the better choice unless there is
  something about SysV that you really need to use. The POSIX
  interface is also less cumbersome. An added bonus (on Linux) is
  that the POSIX mqd_t type representing a queue identifier is
  actually a file descriptor. This allows for the extremely
  convenient use of POSIX queues with select, poll and similar
  facilities.

  Basicly, msgget, msgsnd, msgrcv are System V IPC, while mq_open,
  mq_send, mq_receive are POSIX IPC.

  POSIX IPC is designed after System V IPC. So a lot of old systems
  support only System V IPC while new systems begin to support POSIX
  IPC as well. And, because POSIX IPC can learn the advantages and
  disadvantages from System V IPC, POSIX IPC may be designed and
  implemented better. A notable difference is that
  all the POSIX IPC interfaces are thread safe.

  Comments: SysV interfaces are also thread safe. There would be no
  point to them otherwise
  Well, I didn't confirm it but just read it from
  docs.oracle.com/cd/E19455-01/806-4750/6jdqdfltf/index.html that
  SysV IPC may not thread safe: Unlike the System V IPC interfaces,
  the POSIX IPC interfaces are all multithread safe

# IPC benchmarking
# https://stackoverflow.com/questions/2854004/ipc-speed-and-compare
  Pipe benchmark                     36,539 msg/s
  FIFOs (named pipes) benchmark      26,246 msg/s
  Message Queues benchmark           67,920 msg/s
  Shared Memory benchmark         3,821,893 msg/s  100x faster!
  TCP sockets benchmark              22,483 msg/s
  Unix domain sockets benchmark      40,683 msg/s
  ZeroMQ benchmark                   15,414 msg/s

  https://www.qnx.com/developers/docs/6.5.0SP1.update/..
  com.qnx.doc.neutrino_sys_arch/ipc.html#Shared_memory

  Shared memory offers the highest bandwidth IPC available. Once a
  shared-memory object is created, processes with access to the
  object can use pointers to directly read and write into it

*/

// ----------------------------------------------------------------- //
// Messaging class -> automates deletion of message queue!!
// Wb, Dec03, 2003

#define WB_MSG_PRIV 0666  

   typedef struct __wb_msgbuf {
       long    type;
       time_t  time;
       char    text[128];
       double  val;
   } wbMsgBuf;

   const size_t WBMSGSZ = sizeof(wbMsgBuf) - sizeof(long);

class wbMsg {

  public:

    wbMsg (const char* file, int line,
       const char *istr="Message queue", int key_=0
    ){
       info=istr;

       key = (key_>0) ? key_ : getpid();  
       msgflg = WB_MSG_PRIV;

       if ((msqid = msgget(key, IPC_CREAT | msgflg)) < 0) {
           wblog(FL,"ERR msgget returns %d", -msqid);
           exit(-1);
       }

       if (istr && istr[0]) wblog(file,line,"%NMSG create queue "
          "%s (key=%d (0x%lX), qid=%d)", info.data,key,key,msqid
       );
    };

    ~wbMsg () {
       if (info.data && info.data[0]) wblog (FL,"MSG remove queue "
          "%s (key=%d (0x%lX), qid=%d)", info.data,key,key,msqid);
       msgctl (msqid, IPC_RMID, 0);
    };

    int Rec(const char*,int, wbMsgBuf &rbuf,
       unsigned flag=IPC_NOWAIT 
    );

    int Snd(const char*,int, wbMsgBuf &sbuf,
       const char *istr, double val=0., unsigned type=1
    );

    wbstring info;
    int msqid, msgflg;
    key_t key;

  protected:
  private:
};

int wbMsg::Rec (
    const char* file, int line,
    wbMsgBuf &rbuf,
    unsigned flag  
){

    memset(&rbuf,0,sizeof(wbMsgBuf));

    int i = msgrcv(msqid, &rbuf, WBMSGSZ, 1, flag);
    if (i < 0) {
        if (i!=-1)
        wblog(file,line,
           "ERR can't read message port %d (%d)", msqid, i);
        return -i;
    }

    sprintf(str,"%s", ctime(&rbuf.time));
    str[strlen(str)-1] = 0;

    wblog (file,line,"%NNB! RECEIVED MESSAGE %s (%d; 0x%lX %d)\n\n"
     "   time : %s\n"
     "   text : %s\n"
     "   val  : %g\n"
     "   type : %ld%N",
    info.data ? info.data : "", key, key, msqid,
    str, rbuf.text, rbuf.val, rbuf.type);

    return 0;
}

int wbMsg::Snd (
    const char* file, int line,
    wbMsgBuf &sbuf,  
    const char *istr, double val, 
    unsigned type
){

    int i;

    memset(&sbuf,0,sizeof(wbMsgBuf));

    sbuf.type = type;
    sbuf.time = time((time_t*)NULL);
    strncpy(sbuf.text,istr,127);
    sbuf.val  = val;

    i = msgsnd(msqid, &sbuf, WBMSGSZ, IPC_NOWAIT);
    if (i<0) {
        wblog(FL,"ERR can't send to message port %d (err=%d)", msqid, i);
        return -1;
    }

    sprintf(str,"%s", ctime(&sbuf.time));
    str[strlen(str)-1] = 0;

    wblog (file,line,
     "Message sent ...\n"
     "   type  : %ld\n"
     "   time  : %s\n"
     "   text  : %s\n"
     "   val   : %g",
    sbuf.type, str, sbuf.text, sbuf.val);

    return 0;
}

int WbMSGSnd (
    const char* file, int line,
    unsigned key,
    char *text, double val,
    unsigned type=1
){

    wbMsgBuf sbuf;
    int i, msqid;

    if ((msqid = msgget(key, IPC_CREAT | WB_MSG_PRIV)) < 0) {
        wblog(FL,"ERR can't read message port %d (err=%d)", key, -msqid);
        return -1;
    }

    sbuf.type = type;
    sbuf.time = time((time_t*)NULL);
    strcpy (sbuf.text, text);
    sbuf.val  = val;

    i = msgsnd(msqid, &sbuf, WBMSGSZ, IPC_NOWAIT);
    if (i<0) {
        wblog(FL,"ERR can't send to message port %d (err=%d)", msqid, i);
        return -1;
    }

    sprintf (str, "%s", ctime(&sbuf.time));
    str[strlen(str)-1] = 0;

    wblog(file,line,
     "Message sent (%d; 0x%lX %d) ...\n"
     "   type  : %ld\n"
     "   time  : %s\n"
     "   text  : %s\n"
     "   val   : %g", key, key, msqid,
    sbuf.type, str, sbuf.text, sbuf.val);

    return 0;
}

#endif

