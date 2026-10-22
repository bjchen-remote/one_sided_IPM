#include "mex.h"
#include <algorithm>
#include <cmath>
#include <vector>
#include <limits>

static inline double weno_left(const double v[5], double lineScale, double epsilon) {
    const double safeScale = std::max(lineScale,std::sqrt(std::numeric_limits<double>::min()));
    const double v1=v[0]/safeScale, v2=v[1]/safeScale, v3=v[2]/safeScale;
    const double v4=v[3]/safeScale, v5=v[4]/safeScale;
    const double candidate0=(1.0/3.0)*v1-(7.0/6.0)*v2+(11.0/6.0)*v3;
    const double candidate1=-(1.0/6.0)*v2+(5.0/6.0)*v3+(1.0/3.0)*v4;
    const double candidate2=(1.0/3.0)*v3+(5.0/6.0)*v4-(1.0/6.0)*v5;
    const double beta0=(13.0/12.0)*std::pow(v1-2*v2+v3,2)+(1.0/4.0)*std::pow(v1-4*v2+3*v3,2);
    const double beta1=(13.0/12.0)*std::pow(v2-2*v3+v4,2)+(1.0/4.0)*std::pow(v2-v4,2);
    const double beta2=(13.0/12.0)*std::pow(v3-2*v4+v5,2)+(1.0/4.0)*std::pow(3*v3-4*v4+v5,2);
    const double tau5=std::abs(beta0-beta2);
    const double alpha0=0.1*(1+std::pow(tau5/(beta0+epsilon),2));
    const double alpha1=0.6*(1+std::pow(tau5/(beta1+epsilon),2));
    const double alpha2=0.3*(1+std::pow(tau5/(beta2+epsilon),2));
    return safeScale*(alpha0*candidate0+alpha1*candidate1+alpha2*candidate2)/(alpha0+alpha1+alpha2);
}

static double lagrange_weight(const int nodes[6], int target, int j) {
    double w=1.0;
    for(int k=0;k<6;++k) if(k!=j) w*=double(target-nodes[k])/double(nodes[j]-nodes[k]);
    return w;
}

void mexFunction(int nlhs,mxArray*plhs[],int nrhs,const mxArray*prhs[]) {
    if(nlhs!=1 || nrhs!=5) mexErrMsgIdAndTxt("ipm:FastWenoMexArgs","Expected q,a,J,dx,epsilon.");
    const mwSize nLines=mxGetM(prhs[0]), nNodes=mxGetN(prhs[0]);
    if(nNodes<6 || mxGetM(prhs[1])!=nLines || mxGetN(prhs[1])!=nNodes ||
       mxGetM(prhs[2])!=nLines || mxGetN(prhs[2])!=nNodes)
        mexErrMsgIdAndTxt("ipm:FastWenoMexSize","Inputs q,a,J must have equal 2D size.");
    for(int k=0;k<3;++k) if(!mxIsDouble(prhs[k]) || mxIsComplex(prhs[k]))
        mexErrMsgIdAndTxt("ipm:FastWenoMexType","Only real double arrays are supported.");
    const double *q=mxGetDoubles(prhs[0]), *a=mxGetDoubles(prhs[1]), *J=mxGetDoubles(prhs[2]);
    const double dx=mxGetScalar(prhs[3]), epsilon=mxGetScalar(prhs[4]);
    if(!(dx>0 && epsilon>0)) mexErrMsgIdAndTxt("ipm:FastWenoMexScale","dx and epsilon must be positive.");
    plhs[0]=mxCreateDoubleMatrix(nLines,nNodes,mxREAL);
    double *out=mxGetDoubles(plhs[0]);
    const mwSize nExt=nNodes+6;
    const int loNodes[6]={0,1,2,3,4,5};
    const int hiNodes[6]={-5,-4,-3,-2,-1,0};
    double loW[3][6], hiW[3][6];
    for(int t=0;t<3;++t) for(int j=0;j<6;++j) {
        loW[t][j]=lagrange_weight(loNodes,t-3,j);
        hiW[t][j]=lagrange_weight(hiNodes,t+1,j);
    }
    std::vector<double> qLine(nExt),aLine(nExt),jLine(nExt);
    for(mwSize line=0;line<nLines;++line) {
        for(mwSize node=0;node<nNodes;++node) {
            qLine[node+3]=q[line+node*nLines];
            aLine[node+3]=a[line+node*nLines];
            jLine[node+3]=J[line+node*nLines];
        }
        for(int t=0;t<3;++t) {
            double qLo=0,aLo=0,jLo=0,qHi=0,aHi=0,jHi=0;
            for(int k=0;k<6;++k) {
                const mwSize li=line+k*nLines, hi=line+(nNodes-6+k)*nLines;
                qLo+=loW[t][k]*q[li]; aLo+=loW[t][k]*a[li]; jLo+=loW[t][k]*J[li];
                qHi+=hiW[t][k]*q[hi]; aHi+=hiW[t][k]*a[hi]; jHi+=hiW[t][k]*J[hi];
            }
            qLine[t]=qLo; aLine[t]=aLo; jLine[t]=jLo;
            qLine[nNodes+3+t]=qHi; aLine[nNodes+3+t]=aHi; jLine[nNodes+3+t]=jHi;
        }
        double oldFlux=0;
        for(mwSize face=0;face<=nNodes;++face) {
            double localAlpha=0;
            for(int k=0;k<6;++k) localAlpha=std::max(localAlpha,std::abs(aLine[face+k]/jLine[face+k]));
            double positive[6],negative[6],posScale=0,negScale=0;
            for(int k=0;k<6;++k) {
                const double mapped=jLine[face+k]*qLine[face+k];
                const double physical=aLine[face+k]*qLine[face+k];
                positive[k]=0.5*(physical+localAlpha*mapped);
                negative[k]=0.5*(physical-localAlpha*mapped);
                posScale=std::max(posScale,std::abs(positive[k]));
                negScale=std::max(negScale,std::abs(negative[k]));
            }
            const double scaleFloor=std::sqrt(std::numeric_limits<double>::epsilon())*std::max(posScale,negScale);
            posScale=std::max(posScale,scaleFloor); negScale=std::max(negScale,scaleFloor);
            const double pos[5]={positive[0],positive[1],positive[2],positive[3],positive[4]};
            const double neg[5]={negative[5],negative[4],negative[3],negative[2],negative[1]};
            const double flux=weno_left(pos,posScale,epsilon)+weno_left(neg,negScale,epsilon);
            if(face>0) out[line+(face-1)*nLines]=(flux-oldFlux)/dx/J[line+(face-1)*nLines];
            oldFlux=flux;
        }
    }
}
