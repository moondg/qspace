
" http://tabo.aurealsys.com/code/vimrc.html
" see also default formating /usr/share/vim/vim72/syntax/

  set nocompatible
  set ruler             " show the cursor position in statusline
" set number            " show line numbers
  set showcmd
  set tabstop=4         " just affects tab, but not indentation
  set shiftwidth=1
  set expandtab         " expand all tabs to spaces according to shiftwidth

" set fileencodings=utf-8,latin1
  set fileencodings=latin1
  set encoding=latin1

" filetype on           
  set wrap              
  set textwidth=85

" rewrap paragraph starting from cursor position
  map <C-Q> gq}         " NB! C-W is used to jump between split windows!
                        " NB! C-R is used to redu (after u(ndo))

" set nocindent         " turn off indendentation for C source
" set autoindent        " keep indentation on consecutive lines
  set smartindent       " nosmartindent
  set showmatch         " curser jumbs briefly to corresponding brace
  set matchpairs+=<:>   " match, to be used with % 

" set digraph           
" set nodigraph " by default, no digraph

" set paste / set nopaste - do not indent while pasting
  set pastetoggle=<F9>
                        
  set nobackup
" set pdev=Dell_1110    " default printer
  set pdev=farbdose

" see :help printoptions for more details
" default: %<%f%h%m%=Page %N
  set printheader=%<%f%h%m%=%{strftime('Wb,%b%d,%g')\ }\ \ \ \ \ \ \ \ \ \ \ \ Page\ %N
  " %F - filename including full path
  " strftime('%x')  dd/mm/yy  (see also :help strftime)
  " strftime('%X')  hh:mm:ss

  set printoptions=header:5,wrap:y,top:4pc
  " header:0    do not print header (default: 2)
  " header:5    lines to reserve for header
  " top:2pc     top page margin (default: 5pc)
  " wrap:y      wrap long lines
  " syntax:n    turn off syntax hightlighting
  " portrait:n  print landscape

" set vimfino='100,f1

  if has('gui_running')
     set guioptions-=T  " remove the toolbar
     set lines=40       " 40 lines of text instead of 24, perfect for 1024x768
  endif

  if has("gui_running")
      if has("gui_gtk2")
          set guifont=Courier\ New\ 11
      elseif has("x11")
          set guifont=-*-courier-medium-r-normal-*-*-180-*-*-m-*-*
      else
          set guifont=Courier_New:h11:cDEFAULT
      endif
  endif 


" Switch syntax highlighting on, when the terminal has colors 
" Also switch on highlighting the last used search pattern. 
  if &t_Co > 2 || has("gui_running") 
     syntax on 
   " colorscheme pablo "   for modifications, see ~/.vim/after/syntax/c.vim
   " colorscheme default " for modifications, see ~/.vim/after/syntax/c.vim
     hi preproc ctermfg=Blue "brown
     hi constant ctermfg=Darkred " cterm=underline,bold
     hi comment ctermfg=DarkGreen " DarkGreen
     hi statement ctermfg=Darkblue
     hi type ctermfg=Darkblue

     set hlsearch                  " unset with :nohl or <F8> ...
   " map <F8> :set hlsearch! <CR>  " space after <F8> is essential!!
     map <silent> <F8> :nohl <CR>  " temporarily unset hlsearch
                                   " silent = no output to status line
     hi Search ctermbg=Yellow      " LightGray
  endif 
 
" NB! this requires <C-s> to be deactivated for the (x)term
" e.g. xterm pauses(freezes) otherwise, and continues with <C-q>

" Use CTRL-V to do what CTRL-V used to do
" modified from <C-Q> in mswin.vim -- Wb,Sep13,08
  noremap <C-V> <C-V>

   map <silent> <C-Left>          :call ShowPrevious() <CR>
  imap <silent> <C-Left>    <Esc> :call ShowPrevious() <CR>
   map <silent> <C-Right>         :call ShowNext()     <CR>
  imap <silent> <C-Right>   <Esc> :call ShowNext()     <CR>
" NB! also can use <C-L> and <C-H> (like left/right navigation)

" C-Shift?(Up|Down) does not seem to work !??
" and by now also C-Shift-(Left|Right) -> used same L/R as for next navigation
   map <silent> <C-S-Right>       :call DragRight()    <CR>
  imap <silent> <C-S-Right> <Esc> :call DragRight()    <CR>
   map <silent> <C-S-Left>        :call DragLeft()     <CR>
  imap <silent> <C-S-Left>  <Esc> :call DragLeft()     <CR>

   map <silent> <C-S-L>        :call DragRight()    <CR>
  imap <silent> <C-S-L> <Esc>  :call DragRight()    <CR>
   map <silent> <C-S-H>        :call DragLeft()     <CR>
  imap <silent> <C-S-H> <Esc>  :call DragLeft()     <CR>

"  map <C-S-H>           :echo "Hello world." <CR>

" followed by <nr> gt -> goto specific buffer (I like Info / content)
   map <C-I> :buffers<CR>
  imap <C-I> <Esc> :buffers<CR>

   map <C-s> :w<CR>
  imap <C-s> <Esc> :w<CR>

  command -nargs=0 Q :qa
  map Q :qa<CR>
" map <C-w> :q<CR>

  if has("gui_gtk2")
  " computing nodes only have tiny version of vim installed
    set mouse=c
  " set mouse=n
  endif

  " allow change of cursor focus only when in normal or visual mode
  " (i.e. this allows to select text!)

  " a = accept mouse everywhere
  "   - allows scroll; select text = visual mode!)
  "   - shift left click won't do anything can be used to just
  "     select the window without movint the mouse pointer
  "     block selection: <Ctrl+Shift> plus left mouse

  " c=command line mode
  " i=insert mode only (e.g. allows to select text in command mode)
  " n=normal - e.g. don't change focus while in insert mode (can't select!)
  " v=visual - moves the curser within document

  set matchtime=2 " how many tenth of a second to blink matching brackets (5)
  set scrolloff=1 " minimal number of screen lines to kep above/below cursor

" set list " show special characters such as newline=$, ...

" set incsearch         " vim will search for text as you enter

" set foldmethod=syntax
" set foldnestmax=1


" map wb a<C-R>=strftime("Wb %Y-%m-%d %a %I:%M %p")<CR><Esc>
  map wb a<C-R>=strftime("Wb,%b%d,%y")<CR><Esc>
  imap Wb, <C-R>=strftime("Wb,%b%d,%y")<CR>

" see also default formating /usr/share/vim/vim72/syntax/
  autocmd BufEnter *.m   source $HOME/.vim/matlab.vim
" autocmd BufEnter *.p   source /usr/share/vim/vim71/syntax/perl.vim
  autocmd BufEnter *.tex source $HOME/.vim/tex.vim
  autocmd BufEnter *.[ihc]* source $HOME/.vim/cpp.vim

  autocmd BufEnter *.ih  set syntax=cpp
  autocmd BufEnter *.txt source $HOME/.vim/txt.vim

" abbreviations
  source $HOME/.vim/wordlist.vim

  ab mym andreas.weichselbaum@lmu.de

" .........................................................................
" The 'tabline' option allows you to define your preferred way to tab
" pages labels. This isn't easy, thus an example will be given here.
" modified from :help setting-tabline
" Wb,Aug15,08

" set showtabline=1
  set tabline=%!MyTabLine()

" term=bold,underline cterm=bold,underline
  hi TabLine     ctermfg=black ctermbg=white term=underline cterm=underline
  hi TabLineSel  ctermfg=white ctermbg=black term=underline cterm=underline
  hi TabLineFill ctermfg=white ctermbg=white cterm=underline
" hi TabLine 

  function MyTabLine()
    let k = tabpagenr()
    let n = tabpagenr('$')
    let s = ''

    if n<6
       let rr=range(n)
    else
     " let rr=range(max([0,k-2]),min([n-1,k+2]))
       let rr=range(n)
    endif

    for i in rr
     " select the highlighting
       if i + 1 == k
         let s .= '%#TabLineSel#'
       else
         let s .= '%#TabLine#'
       endif

     " set the tab page number (for mouse clicks)
       let s .= '%' . (i + 1) . 'T'

     " the label is made by MyTabLabel()
       let s .= '[' . (i+1) . '] %{MyTabLabel(' . (i + 1) . ')} '

       if i>=4 && k<=4 && i+1<n
          let s .= '...'
          break
       endif
    endfor

  " after the last tab fill with TabLineFill and reset tab page nr
    let s .= '%#TabLineFill#%T'

  " right-align the label to close the current tab page
    if n > 1
    " let s .= '%=%#TabLine#%999Xclose'
      let s .= '%=%#TabLine#%999XX'
    endif

    return s
  endfunction

  function MyTabLabel(n)
    let buflist = tabpagebuflist(a:n)
    let winnr = tabpagewinnr(a:n)

  " check whether file has been modifed
    let modded = ''
    if getbufvar(a:n,'&modified') !=0
       let modded = ' +'
    endif

    return fnamemodify(bufname(buflist[winnr-1]),':t') . modded
  " see help filename-modifiers
  endfunction

  function! ShowNext()
   " if tabpagenr() == tabpagenr("$"); tabm 0; else
     if tabpagenr() != tabpagenr("$")
        exe "tabn"
     endif
  endfunction

  function! ShowPrevious()
   " if tabpagenr() == 1; exe "tabm ".tabpagenr("$"); else
     if tabpagenr() > 1
        exe "tabp"
     endif
  endfunction

" acknowlegement: http://dotfiles.org/~Tinou/.vimrc -- Wb,Sep23,08
  function! DragRight()
   " if tabpagenr() == tabpagenr("$"); tabm 0; else
     if tabpagenr() != tabpagenr("$")
        exe "tabm ".tabpagenr()
     endif
  endfunction

  function! DragLeft()
   " if tabpagenr() == 1; exe "tabm ".tabpagenr("$"); else
     if tabpagenr() > 1
        exe "tabm ".(tabpagenr()-2)
     endif
  endfunction

" .........................................................................

