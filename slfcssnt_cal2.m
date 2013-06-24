function [ outputE, outputV,outputEf,outputA ] = slfcssnt_cal2(delta,gamma0,gamma1,gamma2,gamma3,gamma4,gamma5,gamma6,indc,nT,nB)
                                       
    N=wksp.Nlayer(indc);

    k=wksp.kMin_sh:wksp.kRes_sh:wksp.kMax_sh;

    const=wksp.e*wksp.d/(2*wksp.dielec);
    alpha=0.05;
    A=zeros(1,N);
    AOld=10*ones(1,N);
    %temporary value for a redundant value Ef

    ntotal=-(nT+nB);
    %E(n,k,phi)/V(n,k,phi)
    n=zeros(1,N);
    nOld=ones(1,N);
    
    %specialized for Bilayer. For other layers, get rid of this part
    %----------------------------------------------------------------------
    
    if nT==0
        n(1)=-nB*0.8;
        n(N)=-nB*0.2;
        
    end
    if nB+nT==0
        n(1)=-nB/2;
        n(N)=-nT/2;
    end
    if nB==nT
        n(1)=-nB/2;
        n(N)=-nT/2;
    end
   

    %----------------------------------------------------------------------
    

    [outputE(:,:),outputV(:,:,:),temp]=DiagH_SC(delta,gamma0,gamma1,gamma2,gamma3,gamma4,gamma5,gamma6,A,indc,0,k);
    EfOld=fermi_level(outputE,indc, ntotal,k);

    EfOld
    Ef=EfOld;
    cnt=0;
    th_cnt=0;
    while 1>0
        cnt=cnt+1;

        Ef=fermi_level_v2(outputE,indc, ntotal,Ef,k);

        for i=1:N-1
            sum1=0;
            sum2=0;
            for j=i+1:N
                sum1=sum1+n(j);
            end
            for j=1:i
                sum2=sum2+n(j);
            end
            A(i+1)=A(i)+const*( (nT-nB)+sum1-sum2);
        end
        if mod(N,2)==0
            U=(A(N/2+1)+A(N/2))/2;
            A=A-U;
        else
            U=A(ceil(N/2));
            A=A-U;
        end
        
%         if sqrt((A-AOld)*(A-AOld)')<0.001
        if sqrt((n-nOld)*(n-nOld)')<10^13;
            th_cnt=th_cnt+1;
            if th_cnt>3
                break;end;
        end
        clc;
        AOld=A;       
        fprintf('nT=%.1f, nB=%.1f',nT/10^16,nB/10^16);
        Ef
        n
        nOld
        -(nB+nT)
        sum(n)
        fprintf('----------------count : %d----------------',cnt);
        nOld=n;
        
        
        [outputE(:,:),outputV(:,:,:),n]=DiagH_SC(delta,gamma0,gamma1,gamma2,gamma3,gamma4,gamma5,gamma6,A,indc,Ef,k);
    

        n=nOld*(1-alpha)+n*alpha;
        A=zeros(1,N);
    end
    outputA=A;
    outputEf=Ef;

end