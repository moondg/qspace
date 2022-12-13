function xytight ()

   h = gca;

   zlm = get(h,'ZLimMode');
   zl  = get(h,'ZLim'   );
   fact= diff(get(h,'YLim'))/diff(zl);

   axis tight

   set(h,'ZLim',zl,'ZLimMode', zlm); % 'DataAspectRatio', [1 1 1/fact]

end

