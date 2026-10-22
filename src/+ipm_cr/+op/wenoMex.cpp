#include "mex.h"
#include <algorithm>
#include <cmath>
#include <vector>

// Optional +ipm_cr unit-speed, unit-metric fifth-order WENO-Z prototype.
// [plus,minus] = wenoPairSmoothMex(q,dx,epsilonFloor,inflowDegree,lo,hi[,mode])
// Each ROW is an independent line; derivative nodes run across COLUMNS.
// minus already denotes -D_flux(q,-1), the right-biased q derivative.
// mode = 0: fixed optimal weights (.1,.6,.3), combined into the equivalent
// fifth-upwind derivatives on a shared seven-point, common-reference stencil.
// mode = 1 (also the six-argument default): smooth RMS-normalized WENO-Z.
//
// Algorithm contract: retain the original +ipm_c candidate reconstructions,
// epsilon floor, boundary-degree/mask semantics, and the deliberate
// common-reference ghost subtractions. Mode 1 retains the WENO-Z weights,
// replacing only the six-point max-abs normalization by
// sqrt(1 + mean(q_i^2)). Its overflow-safe factorization
// uses max solely as an algebraically cancelling rescaling, not as the
// mathematical normalization. This removes normalization winner changes
// and its zero-amplitude homogeneous singularity for fixed boundary masks.
// This optional kernel makes no full-RHS smoothness or temporal-order claim.
// Compile without fast-math: the common-reference subtractions are deliberate.

struct FaceWeights { double left[3], right[3]; };

static bool inflow_flag(const mxArray* a, mwSize line, mwSize lines) {
    if (mxGetNumberOfElements(a)!=lines || mxIsComplex(a) || mxIsSparse(a) ||
        !(mxIsLogical(a) || mxIsDouble(a)))
        mexErrMsgIdAndTxt("ipm_cr:SmoothUnitWenoMask","Masks need one logical/0/1 value per line.");
    if (mxIsLogical(a)) return mxGetLogicals(a)[line]!=0;
    const double value=mxGetDoubles(a)[line];
    if (!(value==0 || value==1))
        mexErrMsgIdAndTxt("ipm_cr:SmoothUnitWenoMask","Numeric masks must contain only zero/one.");
    return value!=0;
}

static void lower_weights(int degree, int targetIndex, double weights[6]) {
    // Rows evaluate at -3,-2,-1 relative to the boundary node 0.
    static const double d1[3][2]={{4,-3},{3,-2},{2,-1}};
    static const double d2[3][3]={{10,-15,6},{6,-8,3},{3,-3,1}};
    static const double d5[3][6]={{56,-210,336,-280,120,-21},
        {21,-70,105,-84,35,-6},{6,-15,20,-15,6,-1}};
    std::fill(weights,weights+6,0.0);
    for (int k=0;k<=degree;++k)
        weights[k]=(degree==1 ? d1[targetIndex][k] :
            (degree==2 ? d2[targetIndex][k] : d5[targetIndex][k]));
}

struct Line {
    std::vector<double> value;
    double lowerReference,upperReference,lowerOffset[3],upperOffset[3];
    mwSize nodes;
    explicit Line(mwSize n):value(n+6),nodes(n) {}
    double difference(mwSize index,double reference) const {
        // Do not round a small ghost offset away by adding its large anchor.
        if (index<3) return (lowerReference-reference)+lowerOffset[index];
        if (index>=nodes+3)
            return (upperReference-reference)+upperOffset[index-nodes-3];
        return value[index]-reference;
    }
};

static void weights(const double v[5],double epsilon,double result[3]) {
    // v are local differences / sqrt(1 + the ORIGINAL six-point mean(q_i^2)).
    const double d12=v[0]-v[1],d23=v[1]-v[2];
    const double d34=v[2]-v[3],d45=v[3]-v[4];
    const double b0=(13.0/12.0)*(d12-d23)*(d12-d23)+
        .25*(d12-3*d23)*(d12-3*d23);
    const double b1=(13.0/12.0)*(d23-d34)*(d23-d34)+
        .25*(d23+d34)*(d23+d34);
    const double b2=(13.0/12.0)*(d34-d45)*(d34-d45)+
        .25*(3*d34-d45)*(3*d34-d45);
    const double tau=std::abs(b0-b2);
    const double r0=tau/(b0+epsilon),r1=tau/(b1+epsilon),r2=tau/(b2+epsilon);
    result[0]=.1*(1+r0*r0);result[1]=.6*(1+r1*r1);result[2]=.3*(1+r2*r2);
    const double sum=result[0]+result[1]+result[2];
    for (int k=0;k<3;++k) result[k]/=sum;
}

static double reconstruct(const double v[5],const double w[3]) {
    const double p0=(1.0/3.0)*v[0]-(7.0/6.0)*v[1]+(11.0/6.0)*v[2];
    const double p1=-(1.0/6.0)*v[1]+(5.0/6.0)*v[2]+(1.0/3.0)*v[3];
    const double p2=(1.0/3.0)*v[2]+(5.0/6.0)*v[3]-(1.0/6.0)*v[4];
    return w[0]*p0+w[1]*p1+w[2]*p2;
}

void mexFunction(int nlhs,mxArray* plhs[],int nrhs,const mxArray* prhs[]) {
    if (nlhs!=2 || !(nrhs==6 || nrhs==7))
        mexErrMsgIdAndTxt("ipm_cr:SmoothUnitWenoArgs",
            "Expected [plus,minus]=wenoPairSmoothMex(q,dx,epsilonFloor,inflowDegree,lo,hi[,reconstructionMode]).");
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]) ||
        mxGetNumberOfDimensions(prhs[0])!=2)
        mexErrMsgIdAndTxt("ipm_cr:SmoothUnitWenoField","q must be a full real double matrix.");
    const mwSize lines=mxGetM(prhs[0]),nodes=mxGetN(prhs[0]);
    if (lines==0 || nodes<6)
        mexErrMsgIdAndTxt("ipm_cr:SmoothUnitWenoSize","At least one line and six nodes are required.");
    for (int k=1;k<=3;++k)
        if (!mxIsDouble(prhs[k]) || mxIsComplex(prhs[k]) || mxIsSparse(prhs[k]) ||
            mxGetNumberOfElements(prhs[k])!=1)
            mexErrMsgIdAndTxt("ipm_cr:SmoothUnitWenoScalars","dx, epsilonFloor and degree must be real double scalars.");
    const double dx=mxGetScalar(prhs[1]),floor=mxGetScalar(prhs[2]),degreeValue=mxGetScalar(prhs[3]);
    if (!(dx>0 && std::isfinite(dx) && floor>0 && std::isfinite(floor)))
        mexErrMsgIdAndTxt("ipm_cr:SmoothUnitWenoScalars","dx and epsilonFloor must be positive and finite.");
    if (!(degreeValue==1 || degreeValue==2 || degreeValue==5))
        mexErrMsgIdAndTxt("ipm_cr:SmoothUnitWenoDegree","Inflow degree must be 1, 2 or 5.");
    int reconstructionMode=1;
    if (nrhs==7) {
        if (!mxIsDouble(prhs[6]) || mxIsComplex(prhs[6]) || mxIsSparse(prhs[6]) ||
            mxGetNumberOfElements(prhs[6])!=1)
            mexErrMsgIdAndTxt("ipm_cr:SmoothUnitWenoMode","reconstructionMode must be a real double scalar, 0 or 1.");
        const double mode=mxGetScalar(prhs[6]);
        if (!(mode==0 || mode==1))
            mexErrMsgIdAndTxt("ipm_cr:SmoothUnitWenoMode","reconstructionMode must be 0 (linear optimal) or 1 (smooth WENO-Z).");
        reconstructionMode=static_cast<int>(mode);
    }
    const int inflowDegree=static_cast<int>(degreeValue);
    const double epsilon=std::max(floor,1.0/(512.0*512.0));
    const double* q=mxGetDoubles(prhs[0]);
    plhs[0]=mxCreateDoubleMatrix(lines,nodes,mxREAL);
    plhs[1]=mxCreateDoubleMatrix(lines,nodes,mxREAL);
    double* plus=mxGetDoubles(plhs[0]);double* minus=mxGetDoubles(plhs[1]);
    Line line(nodes);
    std::vector<FaceWeights> face(reconstructionMode==1 ? nodes+1 : 0);
    for (mwSize row=0;row<lines;++row) {
        const int lowerDegree=inflow_flag(prhs[4],row,lines) ? inflowDegree : 5;
        const int upperDegree=inflow_flag(prhs[5],row,lines) ? inflowDegree : 5;
        for (mwSize j=0;j<nodes;++j) {
            const double value=q[row+j*lines];
            if (!std::isfinite(value))
                mexErrMsgIdAndTxt("ipm_cr:SmoothUnitWenoFinite","q must be finite.");
            line.value[j+3]=value;
        }
        line.lowerReference=line.value[3];line.upperReference=line.value[nodes+2];
        for (int t=0;t<3;++t) {
            double lo[6],hi[6];lower_weights(lowerDegree,t,lo);
            // Upper targets +1,+2,+3 reflect lower targets -1,-2,-3.
            lower_weights(upperDegree,2-t,hi);
            double low=0,high=0;
            for (int k=0;k<=lowerDegree;++k)
                low+=lo[k]*(line.value[3+k]-line.lowerReference);
            for (int k=0;k<=upperDegree;++k)
                high+=hi[upperDegree-k]*(line.value[nodes+2-upperDegree+k]-line.upperReference);
            line.lowerOffset[t]=low;line.upperOffset[t]=high;
            line.value[t]=line.lowerReference+low;
            line.value[nodes+3+t]=line.upperReference+high;
        }
        if (reconstructionMode==0) {
            // Combine the optimal face reconstructions analytically:
            // D+ [-2,15,-60,20,30,-3]/(60*h), offsets -3..2;
            // D- [3,-30,-20,60,-15,2]/(60*h), offsets -2..3.
            // The central difference is exactly zero. Share the remaining
            // six differences, preserving the original ghost references.
            for (mwSize j=0;j<nodes;++j) {
                const double reference=line.value[j+3];
                const double dm3=line.difference(j,reference);
                const double dm2=line.difference(j+1,reference);
                const double dm1=line.difference(j+2,reference);
                const double dp1=line.difference(j+4,reference);
                const double dp2=line.difference(j+5,reference);
                const double dp3=line.difference(j+6,reference);
                // Scale the coefficients before summing, avoiding an
                // unnecessary factor-60 intermediate or 60*h overflow.
                plus[row+j*lines]=(-dm3/30.0+dm2/4.0-dm1+dp1/2.0-dp2/20.0)/dx;
                minus[row+j*lines]=(dm2/20.0-dm1/2.0+dp1-dp2/4.0+dp3/30.0)/dx;
                if (!(std::isfinite(plus[row+j*lines]) && std::isfinite(minus[row+j*lines])))
                    mexErrMsgIdAndTxt("ipm_cr:SmoothUnitWenoFinite","A unit derivative is not representable.");
            }
            continue;
        }
        if (reconstructionMode==1) {
            for (mwSize f=0;f<=nodes;++f) {
                // scale = mag*reducedScale = sqrt(1 + sum(q_i*q_i)/6).
                // Divide in two steps instead of forming scale: even when all
                // q_i are representable, their squares or scale may overflow.
                double mag=1.0;
                for (int k=0;k<6;++k) mag=std::max(mag,std::abs(line.value[f+k]));
                if (!std::isfinite(mag))
                    mexErrMsgIdAndTxt("ipm_cr:SmoothUnitWenoFinite","An extrapolated ghost value is not representable.");
                double squareSum=0.0;
                for (int k=0;k<6;++k) {
                    const double reduced=line.value[f+k]/mag;
                    squareSum+=reduced*reduced;
                }
                const double inverseMag=1.0/mag;
                const double reducedScale=std::sqrt(inverseMag*inverseMag+squareSum/6.0);
                const double anchor=line.value[f+2];double left[5],right[5];
                for (int k=0;k<5;++k) {
                    left[k]=(line.difference(f+k,anchor)/mag)/reducedScale;
                    right[k]=(line.difference(f+5-k,anchor)/mag)/reducedScale;
                }
                weights(left,epsilon,face[f].left);weights(right,epsilon,face[f].right);
            }
        }
        for (mwSize j=0;j<nodes;++j) {
            const double reference=line.value[j+3];double left[2],right[2];
            for (int side=0;side<2;++side) {
                const mwSize f=j+side;double l[5],r[5];
                for (int k=0;k<5;++k) {
                    l[k]=line.difference(f+k,reference);
                    r[k]=line.difference(f+5-k,reference);
                }
                left[side]=reconstruct(l,face[f].left);
                right[side]=reconstruct(r,face[f].right);
            }
            plus[row+j*lines]=(left[1]-left[0])/dx;
            minus[row+j*lines]=(right[1]-right[0])/dx;
            if (!(std::isfinite(plus[row+j*lines]) && std::isfinite(minus[row+j*lines])))
                mexErrMsgIdAndTxt("ipm_cr:SmoothUnitWenoFinite","A unit derivative is not representable.");
        }
    }
}
