function plotQSpectra(A,varargin)
% function plotQSpectra(A,varargin)
% plain wrapper routine => see @QSpace/plotQSpectra.m
% Wb,Mar13,13

  if nargin>1 && ischar(A)
     f=sprintf('%s_%02g.mat',A,varargin{1});
     if ~exist(f,'file')
        error('Wb:ERR','\n   ERR invalid file %s',repHome(f)); end
     q=load(f,'HK'); A=QSpace(q.HK); varargin=varargin(2:end);
  elseif nargin && isnumber(A)
     f=sprintf([getenv('LMA') '/NRG/NRG_%02g.mat'],A);
     if ~exist(f,'file')
        error('Wb:ERR','\n   ERR invalid file %s',repHome(f)); end
     q=load(f,'HK'); A=QSpace(q.HK);
  end

  if isset('f'), o={'mat',untex(repHome(f))}; else o={}; end

  plotQSpectra(QSpace(A),varargin{:},o{:});

end

