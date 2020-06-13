% -------------------------------------------------------------------
% Copyright (C) 2006-2011 Ting Chen, Baba C. Vemuri and Anand Rangarajan
% 
% Authors: Ting Chen
% Date:    12/09/2011
% 
% Contact Information:
%
% Ting Chen, tichen@cise.ufl.edu
% 
% Terms:	  
% 
% The source code (M-files) are provided under the
% terms of the GNU General Public License with an explicit
% clause permitting the execution of the M-files from within
% a MATLAB environment. See the LICENSE file for details.
% ------------------------------------------------------------------

function gv=GradP2(Pts1,Pts2,tag)

n1=max(size(Pts1));
n2=max(size(Pts2));
d=min(size(Pts1));


if tag==1 %Pts1=Pts2
    for i=1:n1
        vx(i)=0;
        vy(i)=0;
    end

for i=1:n1
    for j=1:n2
        if i~=j
            vx(i)=vx(i)+2*min(Pts1(i,2),Pts2(j,2))*gradymin(Pts1(i,1),Pts2(j,1));
            vy(i)=vy(i)+2*min(Pts1(i,1),Pts2(j,1))*gradymin(Pts1(i,2),Pts2(j,2));
        end
        if i==j
            vx(i)=vx(i)+min(Pts1(i,2),Pts2(j,2))*1;
            vy(i)=vy(i)+min(Pts1(i,1),Pts2(j,1))*1;
        end
    end
end

gv=1/(n1*n2)*[vx',vy'];

end
%--------------------------------------------------------------%
if tag==2 %Pts1!=Pts2
    
for i=1:n1
    vx(i)=0;
    vy(i)=0;
end
    
for i=1:n1
    for j=1:n2
        vx(i)=vx(i)+gradymin(Pts1(i,1),Pts2(j,1))*min(Pts1(i,2),Pts2(j,2));
        vy(i)=vy(i)+min(Pts1(i,1),Pts2(j,1))*gradymin(Pts1(i,2),Pts2(j,2));
    end
end

gv=1/(n1*n2)*[vx',vy'];

end
