function strain_gauage
lam=1;% lame constant lambda
miu=1;%lame constant miu
h=1; %height of the left and the right part
l=1; %length of the middle part
w=0.1*h; %width of the part
u0=0.001*h; %displacement
c1=[1,0,0,0;0,1/2,1/2,0;0,0,0,1;0,1/2,1/2,0];
c2=[lam+2*miu,0,lam,0;0,miu,0,miu;lam,0,lam+2*miu,0;0,miu,0,miu];
c=c2*c1;
c11=c(1:2,1:2);
c12=c(1:2,3:4);
c21=c(3:4,1:2);
c22=c(3:4,3:4);
c=[c11(:);c21(:);c12(:);c22(:)]; %this part gives the coefficient of the 2D elastic equation
model=createpde(2);
leftrect=[3 4 0 w w 0 0 0 h h];
uprect=[3 4 0 0 l l h h-w h-w h];
rightrect=[3 4 l-w l-w l l h 0 0 h];
gdm=[leftrect;uprect;rightrect]';
[dl,bt]=decsg(gdm,'R1+R2+R3',['R1';'R2';'R3']');
g=csgdel(dl,bt);
geometryFromEdges(model,g);

figure(1);
pdegplot(model, 'EdgeLabels', 'on', 'FaceLabels', 'on','SubdomainLabels','on'); 
xlabel('X-coordinate, meters')
ylabel('Y-coordinate, meters')
hold on %create the model, and plot the shape of the gauage

applyBoundaryCondition(model,'dirichlet','Edge',[3],...
                       'u',[-u0,0],'EquationIndex',[1,2]);
applyBoundaryCondition(model,'dirichlet','Edge',[4],...
                       'u',[u0,0],'EquationIndex',[1,2]);
applyBoundaryCondition(model,'neumann','Edge',[1,2,5,6,7,8,9,10,11,12]);

specifyCoefficients(model,'m',0,'d',0,'c',c,'a',0,'f',[0;0]);
%boundary conditions and the coefficients
hmax=l/50; %the size of the mesh for the FEM
generateMesh(model,'Hmax',hmax);

result = solvepde(model);
u = result.NodalSolution; %get the result, u(:,1) is the displacement on x direction, u(:,2)is the displacement on z direction
figure(2)
pdeplot(model,'XYData',u,'ZData',u(:,1),'Mesh','on')
xlabel('x')
ylabel('y')
hold off %plot the result