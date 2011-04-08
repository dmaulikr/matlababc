function [SolEt,SolHr,SolHz] = wguide(SolEt,SolHr,SolHz)

SolEt(end,end-1)=SolEt(end-1,end-1);
SolEt(end-1,end)=SolEt(end-1,end-1);
SolEt(1,end-1)=SolEt(2,end-1);
SolEt(2,end)=SolEt(2,end-1);

SolHr(end,end-1)=SolHr(end-1,end-1);
SolHr(end-1,end)=SolHr(end-1,end-1);
SolHr(1,end-1)=SolHr(2,end-1);
SolHr(2,end)=SolHr(2,end-1);

SolHz(end,end-1)=SolHz(end-1,end-1);
SolHz(end-1,end)=SolHz(end-1,end-1);
SolHz(1,end-1)=SolHz(2,end-1);
SolHz(2,end)=SolHz(2,end-1);