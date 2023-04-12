%integer number of triangular lattices
prompt = "What value of n ";
n = input(prompt);

%insert values from minimisation
alphain=0.75;
betain=0.31;

alpha=alphain/(1-alphain);
beta=betain/(1-betain);

%for above value of alpha and beta
phiS0=3.147948273;
phiS1=0.000920421;
phiS2=-0.000541301;
phiA0=-2.009730739;
phiA1=-0.693036202;
phiSA2=-0.004595536;

%lattice structure
a0=4;%lattice vector of the umtwisted cell
a=a0*n; %lattice vector of moire cell 
c=15;
moir = spinw();
moir.genlattice('lat_const', [a a c], 'angled', [90 90 120]);

%the fractional positions of the atoms in the untwisted bilayer
Aa=0;
Ab=0;
Ba=2/3;
Bb=1/3;

%q vectors
    %first harmonic
q1=2*pi/a0 *[-1/3,1/sqrt(3)];
q2=2*pi/a0 *[-1/3,-1/sqrt(3)];
q3=-(q1+q2);
    %second harmonic
q2h1 = 2*q1 + q3;
q2h2 = 2*q3 + q2;
q2h3 = 2*q2 + q1;


%define empty lists to be filled during the loop
    %spin vectors
spinstr=[];
    %interlayer coupling strength 
AJ=[]; 
BJ=[];

theta=2*asin(1/(2*n));
Jh=-1;%nearest neighbour coupling
d=beta*theta^2;%strength of anisotropy related to beta (and lambda)
%d=dd/5;
J=alpha*theta^2;%strength of coupling, related to alpha (and g)


atomnoBl2=0;%initialise number of atoms added

len=0.2*c;%distance between layers
loop=0;


for i = 0:n-1
    %calculate the new fractional positions 
    Aam=Aa/n +i/n;
    Bam=Ba/n +i/n;
    

    for j = 0:n-1

        loop=loop+1;
       
        Abm=Ab/n +j/n;
        Bbm=Bb/n +j/n;

       %generate lattice
        moir.addatom('r', [Aam Abm  0.4], 'S', 1)
        moir.addatom('r', [Bam Bbm  0.4], 'S', 1,'color','DarkCyan')
    
        moir.addatom('r', [Aam Abm  0.6], 'S', 1)
        moir.addatom('r', [Bam Bbm  0.6], 'S', 1,'color','DarkCyan')

        %create interlayer coupling for each atom 
        %it is sufficient to give them different names for the matrices to
        %be saved each loop
         
            %generate coupling at atoms positions
                %generate the cartesian positions 
        Ax=3*a0*Aam-3*a0*Abm/2;
        Ay=sqrt(3)*a0*Abm/2;
        Ar=[Ax,Ay];

        Bx=3*a0*Bam-3*a0*Bbm/2;
        By=sqrt(3)*a0*Bbm/2;
        Br=[Bx,By];
                %eta evaluated at the cartesian position
        Aeta=(cos(dot(q1,Ar))+cos(dot(q2,Ar))+cos(dot(q3,Ar)));
        Beta=(cos(dot(q1,Br))+cos(dot(q2,Br))+cos(dot(q3,Br)));
                %add the value of the interlayer coupling to a list
        %AJ=[AJ,J*Aeta]
        %BJ=[BJ,J*Beta]

        AJ=J*Aeta;
        BJ=J*Beta;

        %create matrices
             %A bonds
        moir.addmatrix('label',['JcA_' num2str(loop)], 'value', AJ);  % +ve is antiferromagnetic
             %B bonds
        moir.addmatrix('label', ['JcB_' num2str(loop)], 'value', BJ);  % +ve is antiferromagnetic
       
        

       %generate spin structure for each atom
        
            %fourier terms calculated to second order, 
            %so calculate eta function for the second order
        Aeta2h=(cos(dot(q2h1,Ar))+cos(dot(q2h2,Ar))+cos(dot(q2h3,Ar)));
        Beta2h=(cos(dot(q2h1,Br))+cos(dot(q2h2,Br))+cos(dot(q2h3,Br)));

            %define angles
       AphiS=phiS0+phiS1*Aeta+phiS1*Aeta2h;
       AphiA=phiA0+phiA1*Aeta+phiA1*Aeta2h;
       
       Aphi1=1/2*(AphiS+AphiA);
       Aphi2=1/2*(AphiS-AphiA);


       BphiS=phiS0+phiS1*Beta+phiS1*Beta2h;
       BphiA=phiA0+phiA1*Beta+phiA1*Beta2h;

       Bphi1=1/2*(BphiS+BphiA);
       Bphi2=1/2*(BphiS-BphiA);

            %define spin vectors
       SA1=[sin(Aphi1);0;cos(Aphi1)];
       SB1=[sin(Bphi1);0;cos(Bphi1)];
       SA2=[sin(Aphi2);0;cos(Aphi2)];
       SB2=[sin(Bphi2);0;cos(Bphi2)];

            %add the spin vector to a list
       spinstr=[spinstr,SA1,SB1,SA2,SB2]; %4 atoms added to unit cell in each iteration
    end

end

%generating bonds
moir.gencoupling('maxDistance', 5, 'forceNoSym', true);
table=moir.table('bond', 1:40);

% Add interlayer coupling 

        %select correct order bond 
        %this is the bond that has the length of distance between layers
bonds=table.idx(table.length==len);%returns all rows of table with correct length
bondno=bonds(1);
moir.table('bond',1:bondno+1);

%loop to add coupling to bond
%note this can't be done in the other loop as the bonds need to be
%generated when all atoms are created 

for i = 0:n-1
    for j = 0:n-1

        %generate atom number
        atomnoA=atomnoBl2+1;
        atomnoB=atomnoA+1;
        atomnoAl2=atomnoB+1;
        atomnoBl2=atomnoAl2+1;
       
        %generate a string for the atom number
        Astr=['atom_' num2str(atomnoA)]
        Bstr=['atom_' num2str(atomnoB)]
        Astrl2=['atom_' num2str(atomnoAl2)]
        Bstrl2=['atom_' num2str(atomnoBl2)]
        
            %add to A bond
        moir.addcoupling('mat', ['JcA_' num2str(atomnoBl2/4)], 'bond',bondno,'atom',{Astr,Astrl2}); 
            %addto B bonds
        moir.addcoupling('mat', ['JcB_' num2str(atomnoBl2/4)], 'bond',bondno,'atom',{Bstr,Bstrl2});

    end

end

% Add in-plane (honeycomb nearest neighbour) coupling
moir.addmatrix('label', 'J1', 'value', Jh);    % +ve is antiferromagnetic
moir.addcoupling('mat', 'J1', 'bond', 1);

%add anisotropy
moir.addmatrix('value',diag([0 0 -d]),'label','d')
moir.addaniso('d')

% add spin structure
moir.genmagstr('mode', 'direct', 'S', spinstr)

moir.optmagsteep('nRun', 100000)

% define the Q vectors between gamma points
q0 = [0; 0; 1]; Qp={};
Qp{1} = ([  0;   0; 0]+q0).*[n; n; 1];
Qp{2} = ([  1/3;   1/3; 0]+q0).*[n; n; 1];
Qp{3} = ([  1/2;   0; 0]+q0).*[n; n; 1];
Qp{4} = ([  0;   0; 0]+q0).*[n; n; 1];
% Calculates and plots the dispersion along a line in reciprocal space

figure;
if n > 3
    chunks = n^3; % Split calculations into chunks to save memory
else
    chunks = 1;
end 

spec = moir.spinwave(Qp, 'hermit', false, 'sortMode', false, 'optmem', chunks); 
% Plots the spin wave spectrum as real / imaginary dipsersion lines
%sw_plotspec(spec,'mode','auto')
% For red lines, use: 'colormap', [255 0 0]
sw_plotspec(spec,'mode','disp','colorbar',false,'lineWidth',0.1,'colormap',[0 0 0])
ylim([0 ceil(max(spec.omega(:)))+0.5])
% Plots the neutron intensity as real / imaginary dipsersion lines

figure;
sw_plotspec(spec,'mode','disp','colorbar',false,'lineWidth',0.1,'colormap',[0 0 0])
ylim([0 abs(real(ceil(max(spec.omega(:)))))/n])

figure;
%sw_plotspec(spec)
sw_plotspec(sw_egrid(spec,'component','Sperp','imagChk',false,'Evect',0:0.005:7.0),'mode','color','dE',0.02); legend('off')

figure;
%sw_plotspec(spec)
sw_plotspec(sw_egrid(spec,'component','Sperp','imagChk',false,'Evect',0:0.005:7.0),'mode','color','dE',0.02); legend('off')
caxis([0 1])

figure;
%sw_plotspec(spec)
sw_plotspec(sw_egrid(spec,'component','Sperp','imagChk',false,'Evect',0:0.005:7.0),'mode','color','dE',0.02); legend('off')
caxis([0 1])
colormap(jet)


