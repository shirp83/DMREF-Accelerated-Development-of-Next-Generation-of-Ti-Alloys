close all;
clear all;

Beta_Volume(10)= 253685;
Beta_Volume(9) = 144841;
Beta_Volume(8) = 242837;
Beta_Volume(7) = 155700;
Beta_Volume(6) = 241896;
Beta_Volume(5) = 117467;
Beta_Volume(4) = 277482;
Beta_Volume(3) = 254643;
Beta_Volume(2) = 141415;
Beta_Volume(1) = 520871-253685;

Beta_Volume=Beta_Volume/2097152;
for kkk=3:10
        
fname=['ebsdalpha' num2str(kkk) '.txt'];
set(gcf,'PaperpositionMode','Auto')
% specify how to align the x-axis in plots
plotx2south

% specify the different phases
%CS = symmetry('m-3m');
%CS1=symmetry('6/mmm','X||a','Y||b*','Z||c*');
CS={...
    symmetry('m-3m','minearl','Beta'),...
    symmetry('6/mmm','X||a','Y||b*','Z||c*','mineral','Alpha')};
SS = symmetry('-1');   % specimen symmetry


% import ebsd data
ebsd = loadEBSD(fname,CS,SS,'interface','generic',...
  'ColumnNames', ...
      { 'id' 'Phase' 'x' 'y'  ...
        'Euler 1' 'Euler 2' 'Euler 3'}...
        ,...
  'Bunge')

  o = get(ebsd('Alpha'),'orientations');
  setMTEXpref('textInterpreter','latex')

 setMTEXpref('textInterpreter','latex')
psi = kernel('de la Vallee Poussin','HALFWIDTH',5*degree)

odf=calcODF(ebsd('Alpha'),'kernel',psi);
figure('position',[100,100,640,480])
plotpdf(odf,Miller(0,0,0,1),'antipodal','silent')
% annotate([xvector,yvector,zvector],'label',{'X','Y','Z'},...
%   'BackgroundColor','w');
annotate([xvector,yvector,zvector],'label',{'X','Y','Z'},'fontname','Times','backgroundcolor','w');
[maxODF,centerODF] = max(odf)
%annotate(centerODF,'label','$q_0$','marker','s','MarkerSize',5,'MarkerFaceColor','r','color','R')
set(gcf,'PaperpositionMode','Auto')
epsname=['alphaodflast_' num2str(kkk) '.eps'];
print(gcf,'-depsc2','-cmyk','-r300','-painters',epsname);

end
