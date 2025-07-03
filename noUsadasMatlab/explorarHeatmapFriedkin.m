% explorarHeatmapFriedkin.m
function [C,P,I] = explorarHeatmapFriedkin(n,maxT,lambVec,nSim,nIter,opR,tolC,tolP,seed)
    p=0.3; simFlag=true; C=zeros(maxT+1,numel(lambVec)); P=C; I=C; idx=1;
    for t=0:maxT
        for j=1:numel(lambVec)
            [W,~,trol,~]=generarModeloFriedkin(p,simFlag,n,floor(t/2),ceil(t/2),opR,seed+idx);
            lam=lambVec(j); lamb=lam*ones(n,1); lamb(trol)=0;
            res=strings(nSim,1);
            for s=1:nSim
                x0=opR(1)+diff(opR)*rand(n,1);
                xf=simularFriedkin(W,x0,lamb,nIter);
                gap=max(xf(:,end))-min(xf(:,end));
                if gap<tolC, res(s)='consenso';
                elseif abs(median(xf(xf<mean(xf(:,end)))')-median(xf(xf>mean(xf(:,end)))'))>tolP, res(s)='polarizacion';
                else res(s)='indefinido'; end
            end
            ct=countcats(categorical(res,{'consenso','polarizacion','indefinido'}));
            C(t+1,j)=100*ct(1)/nSim; P(t+1,j)=100*ct(2)/nSim; I(t+1,j)=100*ct(3)/nSim;
            idx=idx+1;
        end
    end
end
