tic
ndivx = 100;
ndivy = 40;
nbnd = 3;
M=8;
totnode = ndivx*(ndivy-2*ndivy/5)+ndivx*M*(2*M*ndivy/5);
nt =1000;
maxfam = 2000;



	coord=zeros(totnode,2);
	numfam=zeros(totnode,1);
    numfam2=zeros(totnode,1);
	pointfam=zeros(totnode,1);
    pointfam2=zeros(totnode,1);
	pforce=zeros(totnode,2);
	pforceold=zeros(totnode,2);
	bforce=zeros(totnode,2);
	stendens=zeros(totnode,2);
	fncst=zeros(totnode,2); 
	disp=zeros(totnode,2);
	vel=zeros(totnode,2);
	acc=zeros(totnode,2);
	acconehalf=zeros(totnode,2);
	acctwohalf=zeros(totnode,2);
	accthreehalf=zeros(totnode,2);
	velonehalf=zeros(totnode,2);
	veltwohalf=zeros(totnode,2);
	velthreehalf=zeros(totnode,2);
    disponehalf=zeros(totnode,2);
    disptwohalf=zeros(totnode,2);
    dispthreehalf=zeros(totnode,2);
	disponeold=zeros(totnode,2);
	veloneold=zeros(totnode,2);
	acconeold=zeros(totnode,2);
	disponenew=zeros(totnode,2);
	velonenew=zeros(totnode,2);
	acconenew=zeros(totnode,2);
    findpoint=zeros(totnode,1);
    findpoint2=zeros(totnode,1);
    theta=zeros(totnode,1);
	massvec=zeros(totnode,2);
	fail=zeros(totnode,maxfam);
	dmg=zeros(totnode,1);
	nodefam=zeros(totnode*maxfam,1);
    nodefam2=zeros(totnode*maxfam,1);

length = 0.1d0;
width = 0.04d0;
for i=1:totnode
dx(i,1) = length /ndivx;
delta(i,1) = 3.015 * dx(i,1);
thick = dx(i,1);
%delta(i,1): Horizon
end
for i=totnode/2-M*ndivy*(1*M*ndivx/5)+1:totnode/2+M*ndivy*(1*M*ndivx/5)
     dx(i,1) = length /(M*ndivx) ;
     delta(i,1) = 3.015 * dx(i,1);
     thick = dx(i,1);
% dx(i,1) = length / ndivx;
end
dens = 2235.0d0;
emod = 65.0d9;
pratio = 0.2d0;
G0=204.0d0;
for i=1:totnode
kkk=emod/(2.0d0*(1.0d0-pratio));%体积
mmm=emod/(2.0d0*(1.0d0+pratio));%剪切
aal=(kkk-2.0d0*mmm)/2.0d0;
bbl(i,1)=6.0d0*mmm/(pi*delta(i,1)^4.0d0);
ddl(i,1)=2.0d0/(pi*delta(i,1)^3.0d0);
s0(i,1)=sqrt(G0/((6.0d0/pi*mmm+16.0d0/(9.0d0*pi^2.0d0)*(kkk-2.0d0*mmm))*delta(i,1)));
area(i,1) = dx(i,1) * dx(i,1);
vol(i,1) =area(i,1);
radij(i,1)  = dx(i,1)/2.0d0;
m(i,1)=(2*pi*delta(i,1)^3)/3;
bc(i,1) = 9.0d0 * emod / (pi * thick* (delta(i,1)^3));
end
dt = 2.0d-8;
slimt=dt/M;
totime = nt * dt;
ctime = 0.0d0;
idist = 0.0d0;
for i = 1:nt
	enddisp(i,1) = 0.0d0;
	endtime(i,1) = 0.0d0;
end
fac = 0.0d0;

nnum = 0;
cnode = 0;
nlength  = 0.0d0;
dforce1 = 0.0d0;
dforce2 = 0.0d0;
crlength = 0.05d0;
% s0 = 0.04472d0;
newcrlength = 0.0d0;
maxdmg = 0.5d0;
maxxc = crlength / 2.0d0;
% denote interge m


%Initialization of fail flag array
%1 means no failure, 0 means failure of the PD bond
fail=zeros(totnode,maxfam);
fail2=zeros(totnode,maxfam);
% for i = 1:totnode
% 	for  j = 1:maxfam
% 		fail(i,j) = 1;
%         fail2(i,j) = 1;
%     end
% end

%Specification of the locations of material points
%Material points of the internal region
for i=1:totnode
    findpoint2(i,1)=1;
end

for i=1:totnode
    if(i<=(3*ndivy/10)*ndivx)
        c=fix(i/ndivx)+1;
        d=mod(i,ndivx);
        if(d==0)
            d=ndivx;
            c=c-1;
        end
        if (c<=nbnd)
            findpoint2(i,1)=0;
        end
        test(i,1)=c;
        test(i,2)=d;
        coordx=-1.0d0 /2.0d0 * length + (length/(2*ndivx)) + (d - 1) * length/ndivx;
        coordy = -1.0d0 /2.0d0 * width + (length/(2*ndivx)) + (c- 1) * length/ndivx;  
    end

    if(i>(3*ndivy/10)*ndivx)&&(i<=((3*ndivy/10)*ndivx+M*ndivx*(2*M*ndivy/5)))
        c=fix((i-(3*ndivy/10)*ndivx)/(M*ndivx))+1;
        d=mod((i-(3*ndivy/10)*ndivx),(M*ndivx));
        if(d==0)
            d=M*ndivx;
            c=c-1;
        end
        test(i,1)=c;
        test(i,2)=d;
       coordx = -1.0d0 /2.0d0 * length + (length/(2*M*ndivx)) + (d - 1) * length/(M*ndivx);
       coordy=-1.0d0 /2.0d0 *width+ (length/(2*ndivx)) + ((3*ndivy/10) - 1) * length/ndivx+(length/(2*ndivx))+(length/(2*M*ndivx))+(c-1)*length/(M*ndivx);
        end
    
        
     if(i>((3*ndivx/10)*ndivy+M*ndivy*(2*M*ndivx/5)))
        c=fix((i-((3*ndivx/10)*ndivy+M*ndivy*(2*M*ndivx/5)))/ndivx)+1;
        d=mod((i-((3*ndivx/10)*ndivy+M*ndivy*(2*M*ndivx/5))),ndivx);
        if(d==0)
            d=ndivx;
            c=c-1;
        end
        test(i,1)=c;
        test(i,2)=d;
        if (c>(ndivy*3/10-nbnd)) 
            findpoint2(i,1)=0;
        end
        coordx = -1.0d0 /2.0d0 * length + (length/(2*ndivx)) + (d - 1) * length/ndivx;
        coordy=-1.0d0 /2.0d0 * width + (length/(2*ndivx)) + ((3*ndivy/10) - 1) * length/ndivx+(length/(2*ndivx))+(length/(2*M*ndivx))+((2*M*ndivy/5)-1)*length/(M*ndivx)+(length/(2*M*ndivx))+(length/(2*ndivx))+(c-1)*length/ndivx;
     end
    nnum=nnum+1;
    coord(nnum,1) = coordx;
    coord(nnum,2) = coordy;
end








for i = 1:totnode
    if (i==1)  
        pointfam(i,1) = 1;
    else
        pointfam(i,1) = pointfam(i-1,1) + numfam(i-1,1);
    end
    for j = 1:totnode
        idist = sqrt((coord(j,1) - coord(i,1))^2 + (coord(j,2) - coord(i,2))^2);
        if (i~=j) 
            if(idist<=delta(i,1)) 
                numfam(i,1) = numfam(i,1) + 1;
                nodefam(pointfam(i,1)+numfam(i,1)-1,1) = j;
            end
        end
    end
end


for i = 1:totnode
    if (i==1)  
        pointfam2(i,1) = 1;
    else
        pointfam2(i,1) = pointfam2(i-1,1) + numfam2(i-1,1);
    end
    for j = 1:totnode
        idist = sqrt((coord(j,1) - coord(i,1))^2 + (coord(j,2) - coord(i,2))^2);
        if (i~=j) 
            if(idist<=delta(j,1)) 
                numfam2(i,1) = numfam2(i,1) + 1;
                nodefam2(pointfam2(i,1)+numfam2(i,1)-1,1) = j;
            end
        end
    end
end



 for i=1:totnode
     if(coord(i,2)<=(-width/2.0d0))||(coord(i,2)>=(width/2.0d0)) 
            findpoint(i,1)=0;
     elseif(coord(i,2)<=(-width/5.0d0-6.015*length /ndivx))||(coord(i,2)>=(width/5.0d0+6.015*length /ndivx)) 
            findpoint(i,1)=4;
     elseif(coord(i,2)<=(-width/5.0d0))||(coord(i,2)>=width/5.0d0)  
            findpoint(i,1)=3;
     elseif(coord(i,2)<=(-width/5.0d0+3.015*length/ndivx))||(coord(i,2)>=(width/5.0d0-3.015*length/ndivx))  
            findpoint(i,1)=2;
     else 
            findpoint(i,1)=1;
  end
 end 


for i = 1:totnode
    for  j = 1:numfam(i,1)
        cnode = nodefam(pointfam(i,1)+j-1,1);
        if ((coord(cnode,2) > 0.0d0))&&((coord(i,2) < 0.0d0)) 
            if (coord(i,1)<=1.0d-10) 
				fail(i,j) = 0;
            elseif (coord(cnode,1)<=1.0d-10) 
				fail(i,j) = 0;
            end
       elseif ((coord(i,2) > 0.0d0))&&((coord(cnode,2) < 0.0d0)) 
            if(coord(i,1)<=1.0d-10)  
				fail(i,j) = 0;
			elseif(coord(cnode,1)<=1.0d-10) 
				fail(i,j) = 0;
           end
        end       
    end
end


for i = 1:totnode
    for  j = 1:numfam2(i,1)
        cnode = nodefam2(pointfam2(i,1)+j-1,1);
        if ((coord(cnode,2) > 0.0d0))&&((coord(i,2) < 0.0d0)) 
            if (coord(i,1)<=1.0d-10) 
				fail2(i,j) = 0;
            elseif (coord(cnode,1)<=1.0d-10) 
				fail2(i,j) = 0;
            end
       elseif ((coord(i,2) > 0.0d0))&&((coord(cnode,2) < 0.0d0)) 
            if(coord(i,1)<=1.0d-10)  
				fail2(i,j) = 0;
			elseif(coord(cnode,1)<=1.0d-10) 
				fail2(i,j) = 0;
           end
        end       
    end
end




scrmatrixv=zeros(totnode,maxfam);
scrmatrixw=zeros(totnode,maxfam);
scrmatrixvv=zeros(totnode,maxfam);
scrmatrixvw=zeros(totnode,maxfam);
scrmatrixww=zeros(totnode,maxfam);
for i=1:totnode
  for j=1:numfam(i,1)  
  cnode = nodefam(pointfam(i,1)+j-1,1);
   idist = sqrt((coord(cnode,1) - coord(i,1))^2 + (coord(cnode,2) - coord(i,2))^2);
    wij=1/idist;
    scrmatrixv(i,j)=(2/m(i,1))*wij*area(cnode,1)*(coord(cnode,1) - coord(i,1))/2;
    scrmatrixw(i,j)=(2/m(i,1))*wij*area(cnode,1)*(coord(cnode,2) - coord(i,2))/2;
    scrmatrixvv(i,j)=(16*mmm/m(i,1))*wij*area(cnode,1)*(coord(cnode,1) - coord(i,1))^2/(idist)^2/2;
    scrmatrixvw(i,j)=(16*mmm/m(i,1))*wij*area(cnode,1)*(coord(cnode,1)-coord(i,1))*(coord(cnode,2)-coord(i,2))/(idist)^2/2;
    scrmatrixww(i,j)=(16*mmm/m(i,1))*wij*area(cnode,1)*(coord(cnode,2)-coord(i,2))^2/(idist)^2/2;
  end
end
  


for i=1:totnode
  for j=1:numfam2(i,1)  
  cnode = nodefam2(pointfam2(i,1)+j-1,1);
   idist = sqrt((coord(cnode,1) - coord(i,1))^2 + (coord(cnode,2) - coord(i,2))^2);             
    wij=1/idist;
    scrmatrixv2(i,j)=(2/m(cnode,1))*wij*area(cnode,1)*(coord(cnode,1) - coord(i,1))/2;
    scrmatrixw2(i,j)=(2/m(cnode,1))*wij*area(cnode,1)*(coord(cnode,2) - coord(i,2))/2;
    scrmatrixvv2(i,j)=(16*mmm/m(cnode,1))*wij*area(cnode,1)*(coord(cnode,1) - coord(i,1))^2/(idist)^2/2;
    scrmatrixvw2(i,j)=(16*mmm/m(cnode,1))*wij*area(cnode,1)*(coord(cnode,1)-coord(i,1))*(coord(cnode,2)-coord(i,2))/(idist)^2/2;
    scrmatrixww2(i,j)=(16*mmm/m(cnode,1))*wij*area(cnode,1)*(coord(cnode,2)-coord(i,2))^2/(idist)^2/2;
  end
  end

    
%Initialization of displacements and velocities
for i = 1:totnode
    vel(i,1) = 0.0d0;
    vel(i,2) = 0.0d0;
    disp(i,1) = 0.0d0;
    disp(i,2) = 0.0d0;        
end

tic
   for tt = 1:nt
   tt
    for i = 1:totnode
        dmgpar1 = 0.0d0;
        dmgpar2 = 0.0d0;
          if(findpoint2(i,1)==1) 
        for  j = 1:numfam(i,1)            
                cnode = nodefam(pointfam(i,1)+j-1,1);
                idist = sqrt((coord(cnode,1) - coord(i,1))^2 + (coord(cnode,2) - coord(i,2))^2);
                nlength = sqrt((coord(cnode,1) + disp(cnode,1) - coord(i,1) - disp(i,1))^2 + (coord(cnode,2) + disp(cnode,2) - coord(i,2) - disp(i,2))^2);

               if (idist<=delta(i,1)-radij(i,1)) 
                    fac(i,1) = 1.0d0;
                elseif (idist<=delta(i,1)+radij(i,1)) 
                    fac(i,1) = (delta(i,1)+radij(i,1)-idist)/(2.0d0*radij(i,1));
                else
                    fac(i,1) = 0.0d0;
                end

                
                if (abs((nlength - idist) / idist) > s0(i,1)) 
						fail(i,j) = 0 ;
                end					 
                dmgpar1 = dmgpar1 + fail(i,j) * area(i,1)*fac(i,1);
                dmgpar2 = dmgpar2 + area(i,1)*fac(i,1);           
        end
              end
        dmg(i,1) = 1.0d0 - dmgpar1 / dmgpar2;
		if ((dmg(i,1)>maxdmg))&&((coord(i,1)>maxxc)) 
			maxxc = coord(i,1);
		end
    end

    
      for i = 1:totnode
          if(findpoint2(i,1)==1) 
        for  j = 1:numfam2(i,1)            
                cnode = nodefam2(pointfam2(i,1)+j-1,1);
                idist = sqrt((coord(cnode,1) - coord(i,1))^2 + (coord(cnode,2) - coord(i,2))^2);
                nlength = sqrt((coord(cnode,1) + disp(cnode,1) - coord(i,1) - disp(i,1))^2 + (coord(cnode,2) + disp(cnode,2) - coord(i,2) - disp(i,2))^2);

               if (idist<=delta(i,1)-radij(i,1)) 
                    fac(i,1) = 1.0d0;
                elseif (idist<=delta(i,1)+radij(i,1)) 
                    fac(i,1) = (delta(i,1)+radij(i,1)-idist)/(2.0d0*radij(i,1));
                else
                    fac(i,1) = 0.0d0;
                end

                
                if (abs((nlength - idist) / idist) > s0(i,1)) 
						fail2(i,j) = 0 ;
                end					          
        end
          end
    end

 for i = 1:ndivx
   bforce(i,2) = -24.0d0*1.0d6/dx(i,1);
 end

 for i = totnode-ndivx+1:totnode    
   bforce(i,2) = 24.0d0*1.0d6/dx(i,1);
 end

      if(tt==1) 

   for i= 1:totnode
    theta(i,1)=0.0d0;
        for  j = 1:numfam(i,1)            
                cnode = nodefam(pointfam(i,1)+j-1,1);
                if (fail(i,j)==1) 
                theta(i,1) = theta(i,1)+(scrmatrixv(i,j)*(disp(cnode,1) - disp(i,1))+scrmatrixw(i,j)*(disp(cnode,2) - disp(i,2)));
                else
                theta(i,1) = theta(i,1)+0.0d0;
                end
          end
     end

     
  for i= 1:totnode
        for  j = 1:numfam2(i,1)            
                cnode = nodefam2(pointfam2(i,1)+j-1,1);
                if (fail2(i,j)==1) 
                theta(i,1) = theta(i,1)+scrmatrixv2(i,j)*(disp(cnode,1) - disp(i,1))+scrmatrixw2(i,j)*(disp(cnode,2) - disp(i,2));
                else
                theta(i,1) = theta(i,1)+0.0d0;
                end
          end
  end



    for i = 1:totnode
        pforce(i,1) = 0.0d0;
        pforce(i,2) = 0.0d0;
        for  j = 1:numfam(i,1)            
                cnode = nodefam(pointfam(i,1)+j-1,1);
                if (fail(i,j)==1) 
                    dforce1 = (((kkk-2*mmm)/m(i,1))/(2/m(i,1)))*scrmatrixv(i,j)*(theta(i,1)-theta(cnode,1))+(scrmatrixvv(i,j)*(disp(cnode,1)-disp(i,1))+scrmatrixvw(i,j)*(disp(cnode,2)-disp(i,2)));
                    dforce2 = (((kkk-2*mmm)/m(i,1))/(2/m(i,1)))*scrmatrixw(i,j)*(theta(i,1)-theta(cnode,1))+(scrmatrixvw(i,j)*(disp(cnode,1)-disp(i,1))+scrmatrixww(i,j)*(disp(cnode,2)-disp(i,2)));
                else
                    dforce1 = 0.0d0;
                    dforce2 = 0.0d0;
                end
                pforce(i,1) = pforce(i,1) + dforce1 ;            
                pforce(i,2) = pforce(i,2) + dforce2 ;  
        end
    end


  
 for i = 1:totnode
        for  j = 1:numfam2(i,1)            
                cnode = nodefam2(pointfam2(i,1)+j-1,1);
                if (fail2(i,j)==1) 
                    dforce1 = (((kkk-2*mmm)/m(cnode,1))/(2/m(cnode,1)))*scrmatrixv2(i,j)*(theta(i,1)-theta(cnode,1))+scrmatrixvv2(i,j)*(disp(cnode,1)-disp(i,1))+scrmatrixvw2(i,j)*(disp(cnode,2)-disp(i,2));
                    dforce2 = (((kkk-2*mmm)/m(cnode,1))/(2/m(cnode,1)))*scrmatrixw2(i,j)*(theta(i,1)-theta(cnode,1))+scrmatrixvw2(i,j)*(disp(cnode,1)-disp(i,1))+scrmatrixww2(i,j)*(disp(cnode,2)-disp(i,2));
                else
                    dforce1 = 0.0d0;
                    dforce2 = 0.0d0;
                end
                pforce(i,1) = pforce(i,1) + dforce1 ;            
                pforce(i,2) = pforce(i,2) + dforce2 ;                      
        end
 end  
 
    for i=1:totnode
      acconehalf(i,1)=(pforce(i,1)+bforce(i,1))/dens;
      acconehalf(i,2)=(pforce(i,2)+bforce(i,2))/dens;
    end


    for i=1:totnode
      velonehalf(i,1)=vel(i,1);
      velonehalf(i,2)=vel(i,2);
    end

    for i=1:totnode
      disponehalf(i,1)=disp(i,1)+2.0d0/3.0d0*velonehalf(i,1)*dt;
      disponehalf(i,2)=disp(i,2)+2.0d0/3.0d0*velonehalf(i,2)*dt;
    end

 

       for i= 1:totnode
    theta(i,1)=0.0d0;
        for  j = 1:numfam(i,1)            
                cnode = nodefam(pointfam(i,1)+j-1,1);
                if (fail(i,j)==1) 
                theta(i,1) = theta(i,1)+(scrmatrixv(i,j)*(disponehalf(cnode,1) - disponehalf(i,1))+scrmatrixw(i,j)*(disponehalf(cnode,2) - disponehalf(i,2)));
                else
                theta(i,1) = theta(i,1)+0.0d0;
                end
          end
     end

     
  for i= 1:totnode
        for  j = 1:numfam2(i,1)            
                cnode = nodefam2(pointfam2(i,1)+j-1,1);
                if (fail2(i,j)==1) 
                theta(i,1) = theta(i,1)+scrmatrixv2(i,j)*(disponehalf(cnode,1) - disponehalf(i,1))+scrmatrixw2(i,j)*(disponehalf(cnode,2) - disponehalf(i,2));
                else
                theta(i,1) = theta(i,1)+0.0d0;
                end
          end
  end



    for i = 1:totnode
        pforce(i,1) = 0.0d0;
        pforce(i,2) = 0.0d0;
        for  j = 1:numfam(i,1)            
                cnode = nodefam(pointfam(i,1)+j-1,1);
                if (fail(i,j)==1) 
                    dforce1 = (((kkk-2*mmm)/m(i,1))/(2/m(i,1)))*scrmatrixv(i,j)*(theta(i,1)-theta(cnode,1))+(scrmatrixvv(i,j)*(disponehalf(cnode,1)-disponehalf(i,1))+scrmatrixvw(i,j)*(disponehalf(cnode,2)-disponehalf(i,2)));
                    dforce2 = (((kkk-2*mmm)/m(i,1))/(2/m(i,1)))*scrmatrixw(i,j)*(theta(i,1)-theta(cnode,1))+(scrmatrixvw(i,j)*(disponehalf(cnode,1)-disponehalf(i,1))+scrmatrixww(i,j)*(disponehalf(cnode,2)-disponehalf(i,2)));
                else
                    dforce1 = 0.0d0;
                    dforce2 = 0.0d0;
                end
                pforce(i,1) = pforce(i,1) + dforce1 ;            
                pforce(i,2) = pforce(i,2) + dforce2 ;  
        end
    end


  
 for i = 1:totnode
        for  j = 1:numfam2(i,1)            
                cnode = nodefam2(pointfam2(i,1)+j-1,1);
                if (fail2(i,j)==1) 
                    dforce1 = (((kkk-2*mmm)/m(cnode,1))/(2/m(cnode,1)))*scrmatrixv2(i,j)*(theta(i,1)-theta(cnode,1))+scrmatrixvv2(i,j)*(disponehalf(cnode,1)-disponehalf(i,1))+scrmatrixvw2(i,j)*(disponehalf(cnode,2)-disponehalf(i,2));
                    dforce2 = (((kkk-2*mmm)/m(cnode,1))/(2/m(cnode,1)))*scrmatrixw2(i,j)*(theta(i,1)-theta(cnode,1))+scrmatrixvw2(i,j)*(disponehalf(cnode,1)-disponehalf(i,1))+scrmatrixww2(i,j)*(disponehalf(cnode,2)-disponehalf(i,2));
                else
                    dforce1 = 0.0d0;
                    dforce2 = 0.0d0;
                end
                pforce(i,1) = pforce(i,1) + dforce1 ;            
                pforce(i,2) = pforce(i,2) + dforce2 ;                      
        end
 end  
    


    for i=1:totnode
      acctwohalf(i,1)=(pforce(i,1)+bforce(i,1))/dens;
      acctwohalf(i,2)=(pforce(i,2)+bforce(i,2))/dens;
    end

    for i=1:totnode
      veltwohalf(i,1)=vel(i,1)+2.0d0/3.0d0*acconehalf(i,1)*dt;
      veltwohalf(i,2)=vel(i,2)+2.0d0/3.0d0*acconehalf(i,2)*dt;
    end

    for i=1:totnode
      disptwohalf(i,1)=disp(i,1)+2.0d0/3.0d0*veltwohalf(i,1)*dt;
      disptwohalf(i,2)=disp(i,2)+2.0d0/3.0d0*veltwohalf(i,2)*dt;
    end



   
       for i= 1:totnode
    theta(i,1)=0.0d0;
        for  j = 1:numfam(i,1)            
                cnode = nodefam(pointfam(i,1)+j-1,1);
                if (fail(i,j)==1) 
                theta(i,1) = theta(i,1)+(scrmatrixv(i,j)*(disptwohalf(cnode,1) - disptwohalf(i,1))+scrmatrixw(i,j)*(disptwohalf(cnode,2) - disptwohalf(i,2)));
                else
                theta(i,1) = theta(i,1)+0.0d0;
                end
          end
     end

     
  for i= 1:totnode
        for  j = 1:numfam2(i,1)            
                cnode = nodefam2(pointfam2(i,1)+j-1,1);
                if (fail2(i,j)==1) 
                theta(i,1) = theta(i,1)+scrmatrixv2(i,j)*(disptwohalf(cnode,1) - disptwohalf(i,1))+scrmatrixw2(i,j)*(disptwohalf(cnode,2) - disptwohalf(i,2));
                else
                theta(i,1) = theta(i,1)+0.0d0;
                end
          end
  end



    for i = 1:totnode
        pforce(i,1) = 0.0d0;
        pforce(i,2) = 0.0d0;
        for  j = 1:numfam(i,1)            
                cnode = nodefam(pointfam(i,1)+j-1,1);
                if (fail(i,j)==1) 
                    dforce1 = (((kkk-2*mmm)/m(i,1))/(2/m(i,1)))*scrmatrixv(i,j)*(theta(i,1)-theta(cnode,1))+(scrmatrixvv(i,j)*(disptwohalf(cnode,1)-disptwohalf(i,1))+scrmatrixvw(i,j)*(disptwohalf(cnode,2)-disptwohalf(i,2)));
                    dforce2 = (((kkk-2*mmm)/m(i,1))/(2/m(i,1)))*scrmatrixw(i,j)*(theta(i,1)-theta(cnode,1))+(scrmatrixvw(i,j)*(disptwohalf(cnode,1)-disptwohalf(i,1))+scrmatrixww(i,j)*(disptwohalf(cnode,2)-disptwohalf(i,2)));
                else
                    dforce1 = 0.0d0;
                    dforce2 = 0.0d0;
                end
                pforce(i,1) = pforce(i,1) + dforce1 ;            
                pforce(i,2) = pforce(i,2) + dforce2 ;  
        end
    end


  
 for i = 1:totnode
        for  j = 1:numfam2(i,1)            
                cnode = nodefam2(pointfam2(i,1)+j-1,1);
                if (fail2(i,j)==1) 
                    dforce1 = (((kkk-2*mmm)/m(cnode,1))/(2/m(cnode,1)))*scrmatrixv2(i,j)*(theta(i,1)-theta(cnode,1))+scrmatrixvv2(i,j)*(disptwohalf(cnode,1)-disptwohalf(i,1))+scrmatrixvw2(i,j)*(disptwohalf(cnode,2)-disptwohalf(i,2));
                    dforce2 = (((kkk-2*mmm)/m(cnode,1))/(2/m(cnode,1)))*scrmatrixw2(i,j)*(theta(i,1)-theta(cnode,1))+scrmatrixvw2(i,j)*(disptwohalf(cnode,1)-disptwohalf(i,1))+scrmatrixww2(i,j)*(disptwohalf(cnode,2)-disptwohalf(i,2));
                else
                    dforce1 = 0.0d0;
                    dforce2 = 0.0d0;
                end
                pforce(i,1) = pforce(i,1) + dforce1 ;            
                pforce(i,2) = pforce(i,2) + dforce2 ;                      
        end
 end  

   
    for i=1:totnode
      accthreehalf(i,1)=(pforce(i,1)+bforce(i,1))/dens;
      accthreehalf(i,2)=(pforce(i,2)+bforce(i,2))/dens;
    end

    for i=1:totnode
      velthreehalf(i,1)=vel(i,1)+2.0d0/3.0d0*acctwohalf(i,1)*dt;
      velthreehalf(i,2)=vel(i,2)+2.0d0/3.0d0*acctwohalf(i,2)*dt;
    end


  for i=1:totnode
        disp(i,1)=disp(i,1)+dt/4.0d0*(velonehalf(i,1)+3.0d0/2.0d0*veltwohalf(i,1)+3.0d0/2.0d0*velthreehalf(i,1));
        disp(i,2)=disp(i,2)+dt/4.0d0*(velonehalf(i,2)+3.0d0/2.0d0*veltwohalf(i,2)+3.0d0/2.0d0*velthreehalf(i,2));
  end

  for i=1:totnode
      vel(i,1)=vel(i,1)+dt/4.0d0*(acconehalf(i,1)+3.0d0/2.0d0*acctwohalf(i,1)+3.0d0/2.0d0*accthreehalf(i,1));
      vel(i,2)=vel(i,2)+dt/4.0d0*(acconehalf(i,2)+3.0d0/2.0d0*acctwohalf(i,2)+3.0d0/2.0d0*accthreehalf(i,2));
   end




for i = 1:totnode
       disponeold(i,1)=0.0d0;
       disponeold(i,2)=0.0d0;
       veloneold(i,1)=0.0d0;
       veloneold(i,2)=0.0d0;
end


  else

%compute coarse region


    for i= 1:totnode
    theta(i,1)=0.0d0;
    if(findpoint(i,1)==3)||(findpoint(i,1)==4) 
        for  j = 1:numfam(i,1)            
                cnode = nodefam(pointfam(i,1)+j-1,1);
                if (fail(i,j)==1) 
                theta(i,1) = theta(i,1)+(scrmatrixv(i,j)*(disp(cnode,1) - disp(i,1))+scrmatrixw(i,j)*(disp(cnode,2) - disp(i,2)));
                else
                theta(i,1) = theta(i,1)+0.0d0;
                end
        end
    end
     end

     
  for i= 1:totnode
      if(findpoint(i,1)==3)||(findpoint(i,1)==4) 
        for  j = 1:numfam2(i,1)            
                cnode = nodefam2(pointfam2(i,1)+j-1,1);
                if (fail2(i,j)==1) 
                theta(i,1) = theta(i,1)+scrmatrixv2(i,j)*(disp(cnode,1) - disp(i,1))+scrmatrixw2(i,j)*(disp(cnode,2) - disp(i,2));
                else
                theta(i,1) = theta(i,1)+0.0d0;
                end
        end
      end
  end



    for i = 1:totnode
        pforce(i,1) = 0.0d0;
        pforce(i,2) = 0.0d0;
        if(findpoint(i,1)==3)||(findpoint(i,1)==4) 
        for  j = 1:numfam(i,1)            
                cnode = nodefam(pointfam(i,1)+j-1,1);
                if (fail(i,j)==1) 
                    dforce1 = (((kkk-2*mmm)/m(i,1))/(2/m(i,1)))*scrmatrixv(i,j)*(theta(i,1)-theta(cnode,1))+(scrmatrixvv(i,j)*(disp(cnode,1)-disp(i,1))+scrmatrixvw(i,j)*(disp(cnode,2)-disp(i,2)));
                    dforce2 = (((kkk-2*mmm)/m(i,1))/(2/m(i,1)))*scrmatrixw(i,j)*(theta(i,1)-theta(cnode,1))+(scrmatrixvw(i,j)*(disp(cnode,1)-disp(i,1))+scrmatrixww(i,j)*(disp(cnode,2)-disp(i,2)));
                else
                    dforce1 = 0.0d0;
                    dforce2 = 0.0d0;
                end
                pforce(i,1) = pforce(i,1) + dforce1 ;            
                pforce(i,2) = pforce(i,2) + dforce2 ;  
        end
        end
    end


%  
 for i = 1:totnode
     if(findpoint(i,1)==3)||(findpoint(i,1)==4) 
        for  j = 1:numfam2(i,1)            
                cnode = nodefam2(pointfam2(i,1)+j-1,1);
                if (fail2(i,j)==1) 
                    dforce1 = (((kkk-2*mmm)/m(cnode,1))/(2/m(cnode,1)))*scrmatrixv2(i,j)*(theta(i,1)-theta(cnode,1))+scrmatrixvv2(i,j)*(disp(cnode,1)-disp(i,1))+scrmatrixvw2(i,j)*(disp(cnode,2)-disp(i,2));
                    dforce2 = (((kkk-2*mmm)/m(cnode,1))/(2/m(cnode,1)))*scrmatrixw2(i,j)*(theta(i,1)-theta(cnode,1))+scrmatrixvw2(i,j)*(disp(cnode,1)-disp(i,1))+scrmatrixww2(i,j)*(disp(cnode,2)-disp(i,2));
                else
                    dforce1 = 0.0d0;
                    dforce2 = 0.0d0;
                end
                pforce(i,1) = pforce(i,1) + dforce1 ;            
                pforce(i,2) = pforce(i,2) + dforce2 ;                      
        end
     end
 end  

    for i=1:totnode
      acconehalf(i,1)=(pforce(i,1) + bforce(i,1))/dens;
      acconehalf(i,2)=(pforce(i,2) + bforce(i,2))/dens;
    end


  for i = 1:totnode
    velonehalf(i,1)=vel(i,1);
    velonehalf(i,2)=vel(i,2);
    end


    for i=1:totnode
      disponehalf(i,1)=disp(i,1)+2.0d0/3.0d0*velonehalf(i,1)*dt;
      disponehalf(i,2)=disp(i,2)+2.0d0/3.0d0*velonehalf(i,2)*dt;
    end
    
%    for i = 1:ndivx*nbnd*3
%         velonehalf(i,2) = -20.0d0;
%         disponehalf(i,2) = -20.0d0 * (tt+3d0/2) * dt;
%     end
% 
%       for i = totnode-ndivx*nbnd*3+1:totnode
%         velonehalf(i,2) = 20.0d0;
%         disponehalf(i,2) = 20.0d0 * (tt+3d0/2) * dt;
%     end 


 for i= 1:totnode
    theta(i,1)=0.0d0;
    if(findpoint(i,1)==3)||(findpoint(i,1)==4) 
        for  j = 1:numfam(i,1)            
                cnode = nodefam(pointfam(i,1)+j-1,1);
                if (fail(i,j)==1) 
                theta(i,1) = theta(i,1)+(scrmatrixv(i,j)*(disponehalf(cnode,1) - disponehalf(i,1))+scrmatrixw(i,j)*(disponehalf(cnode,2) - disponehalf(i,2)));
                else
                theta(i,1) = theta(i,1)+0.0d0;
                end
        end
    end
     end

     
  for i= 1:totnode
      if(findpoint(i,1)==3)||(findpoint(i,1)==4) 
        for  j = 1:numfam2(i,1)            
                cnode = nodefam2(pointfam2(i,1)+j-1,1);
                if (fail2(i,j)==1) 
                theta(i,1) = theta(i,1)+scrmatrixv2(i,j)*(disponehalf(cnode,1) - disponehalf(i,1))+scrmatrixw2(i,j)*(disponehalf(cnode,2) - disponehalf(i,2));
                else
                theta(i,1) = theta(i,1)+0.0d0;
                end
        end
      end
  end



    for i = 1:totnode
        pforce(i,1) = 0.0d0;
        pforce(i,2) = 0.0d0;
        if(findpoint(i,1)==3)||(findpoint(i,1)==4) 
        for  j = 1:numfam(i,1)            
                cnode = nodefam(pointfam(i,1)+j-1,1);
                if (fail(i,j)==1) 
                    dforce1 = (((kkk-2*mmm)/m(i,1))/(2/m(i,1)))*scrmatrixv(i,j)*(theta(i,1)-theta(cnode,1))+(scrmatrixvv(i,j)*(disponehalf(cnode,1)-disponehalf(i,1))+scrmatrixvw(i,j)*(disponehalf(cnode,2)-disponehalf(i,2)));
                    dforce2 = (((kkk-2*mmm)/m(i,1))/(2/m(i,1)))*scrmatrixw(i,j)*(theta(i,1)-theta(cnode,1))+(scrmatrixvw(i,j)*(disponehalf(cnode,1)-disponehalf(i,1))+scrmatrixww(i,j)*(disponehalf(cnode,2)-disponehalf(i,2)));
                else
                    dforce1 = 0.0d0;
                    dforce2 = 0.0d0;
                end
                pforce(i,1) = pforce(i,1) + dforce1 ;            
                pforce(i,2) = pforce(i,2) + dforce2 ;  
        end
        end
    end


%  
 for i = 1:totnode
     if(findpoint(i,1)==3)||(findpoint(i,1)==4) 
        for  j = 1:numfam2(i,1)            
                cnode = nodefam2(pointfam2(i,1)+j-1,1);
                if (fail2(i,j)==1) 
                    dforce1 = (((kkk-2*mmm)/m(cnode,1))/(2/m(cnode,1)))*scrmatrixv2(i,j)*(theta(i,1)-theta(cnode,1))+scrmatrixvv2(i,j)*(disponehalf(cnode,1)-disponehalf(i,1))+scrmatrixvw2(i,j)*(disponehalf(cnode,2)-disponehalf(i,2));
                    dforce2 = (((kkk-2*mmm)/m(cnode,1))/(2/m(cnode,1)))*scrmatrixw2(i,j)*(theta(i,1)-theta(cnode,1))+scrmatrixvw2(i,j)*(disponehalf(cnode,1)-disponehalf(i,1))+scrmatrixww2(i,j)*(disponehalf(cnode,2)-disponehalf(i,2));
                else
                    dforce1 = 0.0d0;
                    dforce2 = 0.0d0;
                end
                pforce(i,1) = pforce(i,1) + dforce1 ;            
                pforce(i,2) = pforce(i,2) + dforce2 ;                      
        end
     end
 end  
    

  for i = 1:totnode
    if(findpoint(i,1)==3)||(findpoint(i,1)==4) 
      acctwohalf(i,1)=(pforce(i,1) + bforce(i,1))/dens;
      acctwohalf(i,2)=(pforce(i,2) + bforce(i,2))/dens;
    end
    end


  for i = 1:totnode
    if(findpoint(i,1)==3)||(findpoint(i,1)==4) 
    veltwohalf(i,1)=vel(i,1)+2.0d0/3.0d0*dt*acconehalf(i,1);
    veltwohalf(i,2)=vel(i,2)+2.0d0/3.0d0*dt*acconehalf(i,2);
    end
    end
 
  for i = 1:totnode
    if(findpoint(i,1)==3)||(findpoint(i,1)==4) 
    disptwohalf(i,1)=disp(i,1)+2.0d0/3.0d0*dt*veltwohalf(i,1);
    disptwohalf(i,2)=disp(i,2)+2.0d0/3.0d0*dt*veltwohalf(i,2);
    end
    end


  for i = 1:totnode
    if(findpoint(i,1)==2) 
      disptwohalf(i,1)=disp(i,1)+2.0d0/3.0d0*dt*vel(i,1)+4.0d0/9.0d0*dt*(vel(i,1)-veloneold(i,1));
      disptwohalf(i,2)=disp(i,2)+2.0d0/3.0d0*dt*vel(i,2)+4.0d0/9.0d0*dt*(vel(i,2)-veloneold(i,2));
    end
  end
    
  
%      for i = 1:ndivx*nbnd*3
%         veltwohalf(i,2) = -20.0d0;
%         disptwohalf(i,2) = -20.0d0 * (tt+3d0/2) * dt;
%     end
% 
%       for i = totnode-ndivx*nbnd*3+1:totnode
%         veltwohalf(i,2) = 20.0d0;
%         disptwohalf(i,2) = 20.0d0 * (tt+3d0/2) * dt;
%     end 

    for i= 1:totnode
    theta(i,1)=0.0d0;
    if(findpoint(i,1)==3)||(findpoint(i,1)==4) 
        for  j = 1:numfam(i,1)            
                cnode = nodefam(pointfam(i,1)+j-1,1);
                if (fail(i,j)==1) 
                theta(i,1) = theta(i,1)+(scrmatrixv(i,j)*(disptwohalf(cnode,1) - disptwohalf(i,1))+scrmatrixw(i,j)*(disptwohalf(cnode,2) - disptwohalf(i,2)));
                else
                theta(i,1) = theta(i,1)+0.0d0;
                end
        end
    end
     end

     
  for i= 1:totnode
      if(findpoint(i,1)==3)||(findpoint(i,1)==4) 
        for  j = 1:numfam2(i,1)            
                cnode = nodefam2(pointfam2(i,1)+j-1,1);
                if (fail2(i,j)==1) 
                theta(i,1) = theta(i,1)+scrmatrixv2(i,j)*(disptwohalf(cnode,1) - disptwohalf(i,1))+scrmatrixw2(i,j)*(disptwohalf(cnode,2) - disptwohalf(i,2));
                else
                theta(i,1) = theta(i,1)+0.0d0;
                end
        end
      end
  end



    for i = 1:totnode
        pforce(i,1) = 0.0d0;
        pforce(i,2) = 0.0d0;
        if(findpoint(i,1)==3)||(findpoint(i,1)==4) 
        for  j = 1:numfam(i,1)            
                cnode = nodefam(pointfam(i,1)+j-1,1);
                if (fail(i,j)==1) 
                    dforce1 = (((kkk-2*mmm)/m(i,1))/(2/m(i,1)))*scrmatrixv(i,j)*(theta(i,1)-theta(cnode,1))+(scrmatrixvv(i,j)*(disptwohalf(cnode,1)-disptwohalf(i,1))+scrmatrixvw(i,j)*(disptwohalf(cnode,2)-disptwohalf(i,2)));
                    dforce2 = (((kkk-2*mmm)/m(i,1))/(2/m(i,1)))*scrmatrixw(i,j)*(theta(i,1)-theta(cnode,1))+(scrmatrixvw(i,j)*(disptwohalf(cnode,1)-disptwohalf(i,1))+scrmatrixww(i,j)*(disptwohalf(cnode,2)-disptwohalf(i,2)));
                else
                    dforce1 = 0.0d0;
                    dforce2 = 0.0d0;
                end
                pforce(i,1) = pforce(i,1) + dforce1 ;            
                pforce(i,2) = pforce(i,2) + dforce2 ;  
        end
        end
    end


%  
 for i = 1:totnode
     if(findpoint(i,1)==3)||(findpoint(i,1)==4) 
        for  j = 1:numfam2(i,1)            
                cnode = nodefam2(pointfam2(i,1)+j-1,1);
                if (fail2(i,j)==1) 
                    dforce1 = (((kkk-2*mmm)/m(cnode,1))/(2/m(cnode,1)))*scrmatrixv2(i,j)*(theta(i,1)-theta(cnode,1))+scrmatrixvv2(i,j)*(disptwohalf(cnode,1)-disptwohalf(i,1))+scrmatrixvw2(i,j)*(disptwohalf(cnode,2)-disptwohalf(i,2));
                    dforce2 = (((kkk-2*mmm)/m(cnode,1))/(2/m(cnode,1)))*scrmatrixw2(i,j)*(theta(i,1)-theta(cnode,1))+scrmatrixvw2(i,j)*(disptwohalf(cnode,1)-disptwohalf(i,1))+scrmatrixww2(i,j)*(disptwohalf(cnode,2)-disptwohalf(i,2));
                else
                    dforce1 = 0.0d0;
                    dforce2 = 0.0d0;
                end
                pforce(i,1) = pforce(i,1) + dforce1 ;            
                pforce(i,2) = pforce(i,2) + dforce2 ;                      
        end
     end
 end  

  for i = 1:totnode
    if(findpoint(i,1)==3)||(findpoint(i,1)==4) 
      accthreehalf(i,1)=(pforce(i,1) + bforce(i,1))/dens;
      accthreehalf(i,2)=(pforce(i,2) + bforce(i,2))/dens;
    end
    end


  for i = 1:totnode
    if(findpoint(i,1)==3)||(findpoint(i,1)==4) 
    velthreehalf(i,1)=vel(i,1)+2.0d0/3.0d0*dt*acctwohalf(i,1);
    velthreehalf(i,2)=vel(i,2)+2.0d0/3.0d0*dt*acctwohalf(i,2);
    end
    end


  for i = 1:totnode
    if(findpoint(i,1)==3)||(findpoint(i,1)==4) 
       disponenew(i,1)=disp(i,1)+dt/4.0d0*(velonehalf(i,1)+3.0d0/2.0d0*veltwohalf(i,1)+3.0d0/2.0d0*velthreehalf(i,1));
       disponenew(i,2)=disp(i,2)+dt/4.0d0*(velonehalf(i,2)+3.0d0/2.0d0*veltwohalf(i,2)+3.0d0/2.0d0*velthreehalf(i,2));
    end
   end

  for i = 1:totnode
    if(findpoint(i,1)==3)||(findpoint(i,1)==4) 
       velonenew(i,1)=vel(i,1)+dt/4.0d0*(acconehalf(i,1)+3.0d0/2.0d0*acctwohalf(i,1)+3.0d0/2.0d0*accthreehalf(i,1));
       velonenew(i,2)=vel(i,2)+dt/4.0d0*(acconehalf(i,2)+3.0d0/2.0d0*acctwohalf(i,2)+3.0d0/2.0d0*accthreehalf(i,2));
    end
   end



%compute fine region




   for mm=1:M
% chazhi 
   for i=1:totnode
    if(findpoint(i,1)==3) 
        alpha=(disponenew(i,1)-disp(i,1)-dt*vel(i,1))/dt^2.0d0;
        beta=1.0d0/(5.0d0*dt)*(2.0d0*alpha-(vel(i,1)-veloneold(i,1))/dt);
        dispslim(i,mm)=disp(i,1)+(mm-1)*slimt*vel(i,1)+((mm-1)*slimt)^2.0d0*(alpha-dt*beta)+((mm-1)*slimt)^3.0d0*beta;
    end
    end


   for i=1:totnode
    if(findpoint(i,1)==3) 
        alpha=(disponenew(i,2)-disp(i,2)-dt*vel(i,2))/dt^2.0d0;
        beta=1.0d0/(5.0d0*dt)*(2.0d0*alpha-(vel(i,2)-veloneold(i,2))/dt);
        dispslim(i,mm+M)=disp(i,2)+(mm-1)*slimt*vel(i,2)+((mm-1)*slimt)^2.0d0*(alpha-dt*beta)+((mm-1)*slimt)^3.0d0*beta;
    end
    end


    for i= 1:totnode
    theta(i,1)=0.0d0;
    if(findpoint(i,1)==1)||(findpoint(i,1)==2) 
        for  j = 1:numfam(i,1)            
                cnode = nodefam(pointfam(i,1)+j-1,1);
                disp11=disp(i,1);
                disp12=disp(i,2);
               if(findpoint(cnode,1)==3) 
                 disp21=dispslim(cnode,mm);
                 disp22=dispslim(cnode,mm+M);
               else
                 disp21=disp(cnode,1);
                 disp22=disp(cnode,2);
           end
                if (fail(i,j)==1) 
                theta(i,1) = theta(i,1)+(scrmatrixv(i,j)*(disp21- disp11)+scrmatrixw(i,j)*(disp22-disp12));
                else
                theta(i,1) = theta(i,1)+0.0d0;
                end
        end
    end
     end

     
  for i= 1:totnode
      if(findpoint(i,1)==1)||(findpoint(i,1)==2) 
        for  j = 1:numfam2(i,1)            
                cnode = nodefam2(pointfam2(i,1)+j-1,1);
                disp11=disp(i,1);
                disp12=disp(i,2);
               if(findpoint(cnode,1)==3) 
                 disp21=dispslim(cnode,mm);
                 disp22=dispslim(cnode,mm+M);
               else
                 disp21=disp(cnode,1);
                 disp22=disp(cnode,2);
               end
                if (fail2(i,j)==1) 
                theta(i,1) = theta(i,1)+scrmatrixv2(i,j)*(disp21 - disp11)+scrmatrixw2(i,j)*(disp22 - disp12);
                else
                theta(i,1) = theta(i,1)+0.0d0;
                end
        end
      end
  end



    for i = 1:totnode
        pforce(i,1) = 0.0d0;
        pforce(i,2) = 0.0d0;
        if(findpoint(i,1)==1)||(findpoint(i,1)==2) 
        for  j = 1:numfam(i,1)            
               cnode = nodefam(pointfam(i,1)+j-1,1);
                disp11=disp(i,1);
                disp12=disp(i,2);
               if(findpoint(cnode,1)==3) 
                 disp21=dispslim(cnode,mm);
                 disp22=dispslim(cnode,mm+M);
               else
                 disp21=disp(cnode,1);
                 disp22=disp(cnode,2);
               end
                if (fail(i,j)==1) 
                    dforce1 = (((kkk-2*mmm)/m(i,1))/(2/m(i,1)))*scrmatrixv(i,j)*(theta(i,1)-theta(cnode,1))+scrmatrixvv(i,j)*(disp21-disp11)+scrmatrixvw(i,j)*(disp22-disp12);
                    dforce2 = (((kkk-2*mmm)/m(i,1))/(2/m(i,1)))*scrmatrixw(i,j)*(theta(i,1)-theta(cnode,1))+scrmatrixvw(i,j)*(disp21-disp11)+scrmatrixww(i,j)*(disp22-disp12);
                else
                    dforce1 = 0.0d0;
                    dforce2 = 0.0d0;
                end
                pforce(i,1) = pforce(i,1) + dforce1 ;            
                pforce(i,2) = pforce(i,2) + dforce2 ;  
        end
        end
    end


%  
 for i = 1:totnode
     if(findpoint(i,1)==1)||(findpoint(i,1)==2) 
        for  j = 1:numfam2(i,1)            
                cnode = nodefam2(pointfam2(i,1)+j-1,1);
                disp11=disp(i,1);
                disp12=disp(i,2);
               if(findpoint(cnode,1)==3) 
                 disp21=dispslim(cnode,mm);
                 disp22=dispslim(cnode,mm+M);
               else
                 disp21=disp(cnode,1);
                 disp22=disp(cnode,2);
               end
                if (fail2(i,j)==1) 
                    dforce1 = (((kkk-2*mmm)/m(cnode,1))/(2/m(cnode,1)))*scrmatrixv2(i,j)*(theta(i,1)-theta(cnode,1))+scrmatrixvv2(i,j)*(disp21-disp11)+scrmatrixvw2(i,j)*(disp22-disp12);
                    dforce2 = (((kkk-2*mmm)/m(cnode,1))/(2/m(cnode,1)))*scrmatrixw2(i,j)*(theta(i,1)-theta(cnode,1))+scrmatrixvw2(i,j)*(disp21-disp11)+scrmatrixww2(i,j)*(disp22-disp12);
                else
                    dforce1 = 0.0d0;
                    dforce2 = 0.0d0;
                end
                pforce(i,1) = pforce(i,1) + dforce1 ;            
                pforce(i,2) = pforce(i,2) + dforce2 ;                      
        end
     end
 end  
    



    for i = 1:totnode
    if(findpoint(i,1)==1)||(findpoint(i,1)==2) 
      acconehalf(i,1)=(pforce(i,1) + bforce(i,1))/dens;
      acconehalf(i,2)=(pforce(i,2) + bforce(i,2))/dens;
    end
    end


    for i = 1:totnode
    if(findpoint(i,1)==1)||(findpoint(i,1)==2) 
      velonehalf(i,1)=vel(i,1);
      velonehalf(i,2)=vel(i,2);
    end
    end

    for i = 1:totnode
    if(findpoint(i,1)==1)||(findpoint(i,1)==2) 
      disponehalf(i,1)=disp(i,1)+2.0d0/3.0d0*velonehalf(i,1)*slimt;
      disponehalf(i,2)=disp(i,2)+2.0d0/3.0d0*velonehalf(i,2)*slimt;
    end
    end
    
    



   for i=1:totnode
    if(findpoint(i,1)==3) 
        alpha=(disponenew(i,1)-disp(i,1)-dt*vel(i,1))/dt^2.0d0;
        beta=1.0d0/(5.0d0*dt)*(2.0d0*alpha-(vel(i,1)-veloneold(i,1))/dt);
        disponehalf(i,1)=dispslim(i,mm)+2.0d0/3.0d0*slimt*(vel(i,1)+((mm-1)*slimt)*2.0d0*(alpha-dt*beta)+3.0d0*((mm-1)*slimt)^2.0d0*beta);
    end
    end


   for i=1:totnode
    if(findpoint(i,1)==3) 
        alpha=(disponenew(i,2)-disp(i,2)-dt*vel(i,2))/dt^2.0d0;
        beta=1.0d0/(5.0d0*dt)*(2.0d0*alpha-(vel(i,2)-veloneold(i,2))/dt);
        disponehalf(i,2)=dispslim(i,mm+M)+2.0d0/3.0d0*slimt*(vel(i,2)+((mm-1)*slimt)*2.0d0*(alpha-dt*beta)+3.0d0*((mm-1)*slimt)^2.0d0*beta);
    end
    end


     for i= 1:totnode
    theta(i,1)=0.0d0;
    if(findpoint(i,1)==1)||(findpoint(i,1)==2) 
        for  j = 1:numfam(i,1)            
                cnode = nodefam(pointfam(i,1)+j-1,1);
                if (fail(i,j)==1) 
                theta(i,1) = theta(i,1)+(scrmatrixv(i,j)*(disponehalf(cnode,1) - disponehalf(i,1))+scrmatrixw(i,j)*(disponehalf(cnode,2) - disponehalf(i,2)));
                else
                theta(i,1) = theta(i,1)+0.0d0;
                end
        end
    end
     end

     
  for i= 1:totnode
      if(findpoint(i,1)==1)||(findpoint(i,1)==2) 
        for  j = 1:numfam2(i,1)            
                cnode = nodefam2(pointfam2(i,1)+j-1,1);
                if (fail2(i,j)==1) 
                theta(i,1) = theta(i,1)+scrmatrixv2(i,j)*(disponehalf(cnode,1) - disponehalf(i,1))+scrmatrixw2(i,j)*(disponehalf(cnode,2) - disponehalf(i,2));
                else
                theta(i,1) = theta(i,1)+0.0d0;
                end
        end
      end
  end



    for i = 1:totnode
        pforce(i,1) = 0.0d0;
        pforce(i,2) = 0.0d0;
        if(findpoint(i,1)==1)||(findpoint(i,1)==2) 
        for  j = 1:numfam(i,1)            
                cnode = nodefam(pointfam(i,1)+j-1,1);
                if (fail(i,j)==1) 
                    dforce1 = (((kkk-2*mmm)/m(i,1))/(2/m(i,1)))*scrmatrixv(i,j)*(theta(i,1)-theta(cnode,1))+(scrmatrixvv(i,j)*(disponehalf(cnode,1)-disponehalf(i,1))+scrmatrixvw(i,j)*(disponehalf(cnode,2)-disponehalf(i,2)));
                    dforce2 = (((kkk-2*mmm)/m(i,1))/(2/m(i,1)))*scrmatrixw(i,j)*(theta(i,1)-theta(cnode,1))+(scrmatrixvw(i,j)*(disponehalf(cnode,1)-disponehalf(i,1))+scrmatrixww(i,j)*(disponehalf(cnode,2)-disponehalf(i,2)));
                else
                    dforce1 = 0.0d0;
                    dforce2 = 0.0d0;
                end
                pforce(i,1) = pforce(i,1) + dforce1 ;            
                pforce(i,2) = pforce(i,2) + dforce2 ;  
        end
        end
    end


%  
 for i = 1:totnode
     if(findpoint(i,1)==1)||(findpoint(i,1)==2) 
        for  j = 1:numfam2(i,1)            
                cnode = nodefam2(pointfam2(i,1)+j-1,1);
                if (fail2(i,j)==1) 
                    dforce1 = (((kkk-2*mmm)/m(cnode,1))/(2/m(cnode,1)))*scrmatrixv2(i,j)*(theta(i,1)-theta(cnode,1))+scrmatrixvv2(i,j)*(disponehalf(cnode,1)-disponehalf(i,1))+scrmatrixvw2(i,j)*(disponehalf(cnode,2)-disponehalf(i,2));
                    dforce2 = (((kkk-2*mmm)/m(cnode,1))/(2/m(cnode,1)))*scrmatrixw2(i,j)*(theta(i,1)-theta(cnode,1))+scrmatrixvw2(i,j)*(disponehalf(cnode,1)-disponehalf(i,1))+scrmatrixww2(i,j)*(disponehalf(cnode,2)-disponehalf(i,2));
                else
                    dforce1 = 0.0d0;
                    dforce2 = 0.0d0;
                end
                pforce(i,1) = pforce(i,1) + dforce1 ;            
                pforce(i,2) = pforce(i,2) + dforce2 ;                      
        end
     end
 end  



    for i = 1:totnode
    if(findpoint(i,1)==1)||(findpoint(i,1)==2) 
      acctwohalf(i,1)=(pforce(i,1) + bforce(i,1))/dens;
      acctwohalf(i,2)=(pforce(i,2) + bforce(i,2))/dens;
    end
    end

    for i = 1:totnode
    if(findpoint(i,1)==1)||(findpoint(i,1)==2) 
      veltwohalf(i,1)=vel(i,1)+2.0d0/3.0d0*acconehalf(i,1)*slimt;
      veltwohalf(i,2)=vel(i,2)+2.0d0/3.0d0*acconehalf(i,2)*slimt;
    end
    end

    for i = 1:totnode
    if(findpoint(i,1)==1)||(findpoint(i,1)==2) 
      disptwohalf(i,1)=disp(i,1)+2.0d0/3.0d0*veltwohalf(i,1)*slimt;
      disptwohalf(i,2)=disp(i,2)+2.0d0/3.0d0*veltwohalf(i,2)*slimt;
    end
    end
    
    



   for i=1:totnode
    if(findpoint(i,1)==3) 
        alpha=(disponenew(i,1)-disp(i,1)-dt*vel(i,1))/dt^2.0d0;
        beta=1.0d0/(5.0d0*dt)*(2.0d0*alpha-(vel(i,1)-veloneold(i,1))/dt);
        disptwohalf(i,1)=disponehalf(i,1)+4.0d0/9.0d0*slimt^2.0d0*(2.0d0*(alpha-dt*beta)+6.0d0*(mm-1)*slimt*beta);
    end
    end


   for i=1:totnode
    if(findpoint(i,1)==3) 
        alpha=(disponenew(i,2)-disp(i,2)-dt*vel(i,2))/dt^2.0d0;
        beta=1.0d0/(5.0d0*dt)*(2.0d0*alpha-(vel(i,2)-veloneold(i,2))/dt);
        disptwohalf(i,2)=disponehalf(i,2)+4.0d0/9.0d0*slimt^2.0d0*(2.0d0*(alpha-dt*beta)+6.0d0*(mm-1)*slimt*beta);
    end
    end



      for i= 1:totnode
    theta(i,1)=0.0d0;
    if(findpoint(i,1)==1)||(findpoint(i,1)==2) 
        for  j = 1:numfam(i,1)            
                cnode = nodefam(pointfam(i,1)+j-1,1);
                if (fail(i,j)==1) 
                theta(i,1) = theta(i,1)+(scrmatrixv(i,j)*(disptwohalf(cnode,1) - disptwohalf(i,1))+scrmatrixw(i,j)*(disptwohalf(cnode,2) - disptwohalf(i,2)));
                else
                theta(i,1) = theta(i,1)+0.0d0;
                end
        end
    end
     end

     
  for i= 1:totnode
      if(findpoint(i,1)==1)||(findpoint(i,1)==2) 
        for  j = 1:numfam2(i,1)            
                cnode = nodefam2(pointfam2(i,1)+j-1,1);
                if (fail2(i,j)==1) 
                theta(i,1) = theta(i,1)+scrmatrixv2(i,j)*(disptwohalf(cnode,1) - disptwohalf(i,1))+scrmatrixw2(i,j)*(disptwohalf(cnode,2) - disptwohalf(i,2));
                else
                theta(i,1) = theta(i,1)+0.0d0;
                end
        end
      end
  end



    for i = 1:totnode
        pforce(i,1) = 0.0d0;
        pforce(i,2) = 0.0d0;
        if(findpoint(i,1)==1)||(findpoint(i,1)==2) 
        for  j = 1:numfam(i,1)            
                cnode = nodefam(pointfam(i,1)+j-1,1);
                if (fail(i,j)==1) 
                    dforce1 = (((kkk-2*mmm)/m(i,1))/(2/m(i,1)))*scrmatrixv(i,j)*(theta(i,1)-theta(cnode,1))+(scrmatrixvv(i,j)*(disptwohalf(cnode,1)-disptwohalf(i,1))+scrmatrixvw(i,j)*(disptwohalf(cnode,2)-disptwohalf(i,2)));
                    dforce2 = (((kkk-2*mmm)/m(i,1))/(2/m(i,1)))*scrmatrixw(i,j)*(theta(i,1)-theta(cnode,1))+(scrmatrixvw(i,j)*(disptwohalf(cnode,1)-disptwohalf(i,1))+scrmatrixww(i,j)*(disptwohalf(cnode,2)-disptwohalf(i,2)));
                else
                    dforce1 = 0.0d0;
                    dforce2 = 0.0d0;
                end
                pforce(i,1) = pforce(i,1) + dforce1 ;            
                pforce(i,2) = pforce(i,2) + dforce2 ;  
        end
        end
    end


%  
 for i = 1:totnode
     if(findpoint(i,1)==1)||(findpoint(i,1)==2) 
        for  j = 1:numfam2(i,1)            
                cnode = nodefam2(pointfam2(i,1)+j-1,1);
                if (fail2(i,j)==1) 
                    dforce1 = (((kkk-2*mmm)/m(cnode,1))/(2/m(cnode,1)))*scrmatrixv2(i,j)*(theta(i,1)-theta(cnode,1))+scrmatrixvv2(i,j)*(disptwohalf(cnode,1)-disptwohalf(i,1))+scrmatrixvw2(i,j)*(disptwohalf(cnode,2)-disptwohalf(i,2));
                    dforce2 = (((kkk-2*mmm)/m(cnode,1))/(2/m(cnode,1)))*scrmatrixw2(i,j)*(theta(i,1)-theta(cnode,1))+scrmatrixvw2(i,j)*(disptwohalf(cnode,1)-disptwohalf(i,1))+scrmatrixww2(i,j)*(disptwohalf(cnode,2)-disptwohalf(i,2));
                else
                    dforce1 = 0.0d0;
                    dforce2 = 0.0d0;
                end
                pforce(i,1) = pforce(i,1) + dforce1 ;            
                pforce(i,2) = pforce(i,2) + dforce2 ;                      
        end
     end
 end  

    for i = 1:totnode
    if(findpoint(i,1)==1)||(findpoint(i,1)==2) 
      accthreehalf(i,1)=(pforce(i,1) + bforce(i,1))/dens;
      accthreehalf(i,2)=(pforce(i,2) + bforce(i,2))/dens;
    end
    end

    for i = 1:totnode
    if(findpoint(i,1)==1)||(findpoint(i,1)==2) 
      velthreehalf(i,1)=vel(i,1)+2.0d0/3.0d0*acctwohalf(i,1)*slimt;
      velthreehalf(i,2)=vel(i,2)+2.0d0/3.0d0*acctwohalf(i,2)*slimt;
    end
    end


 if(mm==1)  
    for i = 1:totnode
    if(findpoint(i,1)==1)||(findpoint(i,1)==2) 
     disponeold(i,1)=disp(i,1);
     disponeold(i,2)=disp(i,2);
     veloneold(i,1)=vel(i,1);
     veloneold(i,2)=vel(i,2);
    end
   end
end


    for i = 1:totnode
    if(findpoint(i,1)==1)||(findpoint(i,1)==2) 
       disp(i,1)=disp(i,1)+1.0d0/4.0d0*slimt*(velonehalf(i,1)+3.0d0/2.0d0*veltwohalf(i,1)+3.0d0/2.0d0*velthreehalf(i,1));
       disp(i,2)=disp(i,2)+1.0d0/4.0d0*slimt*(velonehalf(i,2)+3.0d0/2.0d0*veltwohalf(i,2)+3.0d0/2.0d0*velthreehalf(i,2));
    end
    end


    for i = 1:totnode
    if(findpoint(i,1)==1)||(findpoint(i,1)==2) 
       vel(i,1)=vel(i,1)+1.0d0/4.0d0*slimt*(acconehalf(i,1)+3.0d0/2.0d0*acctwohalf(i,1)+3.0d0/2.0d0*accthreehalf(i,1));
       vel(i,2)=vel(i,2)+1.0d0/4.0d0*slimt*(acconehalf(i,2)+3.0d0/2.0d0*acctwohalf(i,2)+3.0d0/2.0d0*accthreehalf(i,2));
    end
    end



end
 
    for i = 1:totnode
    if(findpoint(i,1)==3)||(findpoint(i,1)==4) 
     disponeold(i,1)=disp(i,1);
     disponeold(i,2)=disp(i,2);
     veloneold(i,1)=vel(i,1);
     veloneold(i,2)=vel(i,2);
    end
   end


   for i=1:totnode
    if(findpoint(i,1)==3)||(findpoint(i,1)==4) 
     disp(i,1)=disponenew(i,1);
     disp(i,2)=disponenew(i,2);
     vel(i,1)=velonenew(i,1);
     vel(i,2)=velonenew(i,2);
    end
    end


end
    end
 
 
toc
dmgtt=zeros(nt,1);
pforce=zeros(nt,1);
    dmgp = reshape(dmg(totnode/2-M*ndivy*(1*M*ndivx/5)+1:totnode/2+M*ndivy*(1*M*ndivx/5)),M*ndivx,2*M*ndivy/5)';

 for i=1:M*ndivx
   X(i)=i;
end

for i=1:2*M*ndivy/5
    Y(i)=i;
end
 
 gca=pcolor(X,Y,dmgp);
 set(gca, 'LineStyle','none');
 colorbar
