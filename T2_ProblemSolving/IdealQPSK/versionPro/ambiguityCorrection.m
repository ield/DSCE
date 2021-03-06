function [rxstream] = ambiguityCorrection(zI, zQ, synchroAmbiguity),
%function [rxstream] = ambiguityCorrection(zI, zQ, synchroAmbiguity),
%Corrects the ambiguity after the pase synchronization
%Input parameters
% zI, zQ: Output bits of the phase synchronizers
%Output parameters
% rxstream: Corrected bits 
%

t=2*synchroAmbiguity(1:2:end)-1; %real part of synchroAmbiguity
signalXcorr_tzI= xcorr(t, zI);
[maxAbsXcorr_tzI, indexI] = max( abs( signalXcorr_tzI  ) );
Xcorr_tzI = signalXcorr_tzI(indexI);

signalXcorr_tzQ= xcorr(t, zQ);
[maxAbsXcorr_tzQ, indexQ] = max( abs( signalXcorr_tzQ  ) );
Xcorr_tzQ = signalXcorr_tzQ(indexQ);    
% Este valor no es maxAbsXcorr porque podr�a ser que el maximo fuera 
% positivo o negativo en function de si el pico esta para arriba o para
% abajo. Dependiendo de hacia donde este ele pico habria un desfase u otro.
% Por ejemplo, si el pico esta positivo en la parte real es como si el
% desfase no hubiera sacado la se�al del cuadrante.


rxstreamI = zeros(size(zI));
rxstreamQ = zeros(size(zQ));
if ( abs( maxAbsXcorr_tzI ) >= abs( maxAbsXcorr_tzQ ) ), 
  if ( Xcorr_tzI>0 ),
     disp('No change')
      rxstreamI(zI<0)=0;
      rxstreamI(zI>=0)=1;
      rxstreamQ(zQ<0)=0;
      rxstreamQ(zQ>=0)=1;
  else
     disp('+180')
     rxstreamI(zI<0)=1;
     rxstreamI(zI>0)=0;
     rxstreamQ(zQ<0)=1;
     rxstreamQ(zQ>0)=0;
  end;
else % abs( maxXcorr_tzI ) < abs( maxXcorr_tzQ )
  if ( Xcorr_tzQ>=0 )
    disp('+90')
    for i=1:length(zI)
      if     (zI(i)>0 && zQ(i)>0), rxstreamI(i) = 1 ; rxstreamQ(i) = 0;
      elseif (zI(i)>0 && zQ(i)<0), rxstreamI(i) = 0 ; rxstreamQ(i) = 0;
      elseif (zI(i)<0 && zQ(i)>0), rxstreamI(i) = 1 ; rxstreamQ(i) = 1;
      else                         rxstreamI(i) = 0 ; rxstreamQ(i) = 1;
      end; %if
    end; %for
  else
    disp('-90')
    for i=1:length(zI)
      if     (zI(i)>0 && zQ(i)>0), rxstreamI(i) = 0 ; rxstreamQ(i) = 1;
      elseif (zI(i)>0 && zQ(i)<0), rxstreamI(i) = 1 ; rxstreamQ(i) = 1;
      elseif (zI(i)<0 && zQ(i)>0), rxstreamI(i) = 0 ; rxstreamQ(i) = 0;
      else                         rxstreamI(i) = 1 ; rxstreamQ(i) = 0;
      end; %if
    end; %for     
  end; %if (maxXcorr_tzQ2
end; %main if

rxstream = zeros(1, 2*length(zI) );
rxstream(1:2:end) = rxstreamI;
rxstream(2:2:end) = rxstreamQ;