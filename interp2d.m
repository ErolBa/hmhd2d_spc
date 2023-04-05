function [out]=interp2d(in,x,z,dir)

% dir=1: smooth along x
%     2: smooth along z
%     3: average both
%
% The output of this can be further smooth by smoothn

[nx nz]=size(in);

% interpolate along x

if (dir==1 || dir==3)

    out1=in;
    
    out1(1,:)=interp1d(in(1,:),z);
    out1(nx,:)=interp1d(in(nx,:),z);


    for k=1:nz
        i=2;
        while in(i,k)==in(1,k) && i<nx
            out1(i,k)=out1(1,k);
            i=i+1;
        end
        i=nx-1;
        while in(i,k)==in(nx,k) && i>1
            out1(i,k)=out1(nx,k);
            i=i-1;
        end
        out1(:,k)=interp1d(out1(:,k),x);
    end
end

if dir==1
    out=out1;
    return
end
% interpolate along z

if (dir==2 || dir==3)
    out2=in;
    
    out2(:,1)=interp1d(in(:,1),x);
    out2(:,nz)=interp1d(in(:,nz),x);


    for i=1:nx
        k=2;
        while in(i,k)==in(i,1) && k<nz
            out2(i,k)=out2(i,1);
            k=k+1;
        end
        k=nz-1;
        while in(i,k)==in(i,nz) && k>1
            out2(i,k)=out2(i,nz);
            k=k-1;
        end
        out2(i,:)=interp1d(out2(i,:),z);
    end

end

if dir==2
    out=out2;
    return
end

out=(out1+out2)/2;
return

function out=interp1d(in,z)

in=in(:);


nz=length(in);
out=zeros(nz,1);
out(1)=in(1);
out(nz)=in(nz);
k=2;
lv=in(1);
lz=z(1);
while k<nz

    k1=k+1;
    while (k1<=nz && in(k1)==in(k))
        k1=k1+1;
    end
    if k1<=nz
        rz=(z(k1-1)+z(k1))/2;
        rv=(in(k1-1)+in(k1))/2;
    else
        rz=z(nz);
        rv=out(nz);
    end
    if (in(k)-lv)*(in(k)-rv)<=0
        out(k:k1-1)=lv+(z(k:k1-1)-lz)*(rv-lv)/(rz-lz);
    else  % extrema
        for k2=k:k1-1
            if z(k2)<(rz+lz)/2
                out(k2)=lv+2*(z(k2)-lz)/(rz-lz)*(in(k)-lv);
            else
                out(k2)=rv+2*(rz-z(k2))/(rz-lz)*(in(k)-rv);
            end
        end
    end

    k=k1;
    lv=rv;
    lz=rz;
end


return

