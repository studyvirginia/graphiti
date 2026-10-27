#pragma once
#include "point_segment_kernel.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>
namespace Sketch {
struct RobustContinuousConfig { double huberDelta{3.0}; double sigmaFloor{1e-6}; double tangentWeight{0.05}; double curvatureWeight{0.01}; };
struct RobustContinuousResult { double fitness{std::numeric_limits<double>::max()}; double position{0}; double tangent{0}; double curvature{0}; bool ok{false}; };
inline auto robustHuber(double x, double d)->double { double a=std::abs(x); return a<=d?0.5*a*a:d*(a-0.5*d); }
inline auto robustSampleLinear(std::vector<Vec2> const& p,double u)->Vec2 { if(p.size()==1)return p[0]; double q=std::clamp(u,0.0,1.0)*(p.size()-1), f=q-std::floor(q); auto i=static_cast<std::size_t>(q),j=std::min(i+1,p.size()-1); return {static_cast<Operon::Scalar>((1-f)*p[i].x+f*p[j].x),static_cast<Operon::Scalar>((1-f)*p[i].y+f*p[j].y)}; }
inline auto robustArcCoordinates(std::vector<Vec2> const& p)->std::vector<double>{
 std::vector<double>a(p.size(),0.0);for(std::size_t i=1;i<p.size();++i)a[i]=a[i-1]+std::hypot(double(p[i].x-p[i-1].x),double(p[i].y-p[i-1].y));
 if(a.empty()||!(a.back()>0.0)||!std::isfinite(a.back()))return{};for(double&v:a)v/=a.back();return a;
}
inline auto robustSampleArc(std::vector<Vec2> const&p,std::vector<double> const&a,double u)->Vec2{
 if(p.size()==1)return p[0];double q=std::clamp(u,0.0,1.0);auto hi=std::lower_bound(a.begin(),a.end(),q);if(hi==a.begin())return p.front();if(hi==a.end())return p.back();auto j=static_cast<std::size_t>(hi-a.begin()),i=j-1;double span=a[j]-a[i],f=span>0.0?(q-a[i])/span:0.0;return{static_cast<Operon::Scalar>((1-f)*p[i].x+f*p[j].x),static_cast<Operon::Scalar>((1-f)*p[i].y+f*p[j].y)};
}
inline auto evaluateRobustContinuous(std::vector<Vec2> const& r,std::vector<Vec2> const& c,double sigma,RobustContinuousConfig const& cfg={})->RobustContinuousResult {
 RobustContinuousResult o; if(r.size()<2||c.size()<2||!std::isfinite(sigma))return o; auto fin=[](Vec2 const&p){return std::isfinite(p.x)&&std::isfinite(p.y);}; for(auto&p:r)if(!fin(p))return o; for(auto&p:c)if(!fin(p))return o;
 auto const rArc=robustArcCoordinates(r),cArc=robustArcCoordinates(c);if(rArc.empty()||cArc.empty())return o;
 double s=std::max(std::abs(sigma),cfg.sigmaFloor),d=std::max(cfg.huberDelta,1e-9); std::vector<Vec2>rr(r.size()),a(r.size());
 double const invH=static_cast<double>(r.size()-1), invH2=invH*invH;
 for(std::size_t i=0;i<r.size();++i){double u=double(i)/double(r.size()-1);rr[i]=robustSampleArc(r,rArc,u);a[i]=robustSampleArc(c,cArc,u); o.position+=robustHuber(std::hypot(double(a[i].x-rr[i].x),double(a[i].y-rr[i].y))/s,d);} o.position/=r.size();
 for(std::size_t i=1;i<r.size();++i)o.tangent+=robustHuber(invH*std::hypot(double((a[i].x-a[i-1].x)-(rr[i].x-rr[i-1].x)),double((a[i].y-a[i-1].y)-(rr[i].y-rr[i-1].y)))/s,d); o.tangent/=r.size()-1;
 if(r.size()>2){for(std::size_t i=1;i+1<r.size();++i)o.curvature+=robustHuber(invH2*std::hypot(double((a[i+1].x-2*a[i].x+a[i-1].x)-(rr[i+1].x-2*rr[i].x+rr[i-1].x)),double((a[i+1].y-2*a[i].y+a[i-1].y)-(rr[i+1].y-2*rr[i].y+rr[i-1].y)))/s,d);o.curvature/=r.size()-2;}
 o.fitness=o.position+cfg.tangentWeight*o.tangent+cfg.curvatureWeight*o.curvature;o.ok=std::isfinite(o.fitness);if(!o.ok)o.fitness=std::numeric_limits<double>::max();return o;
}
inline auto coarseFourPointForwardLowerBound(std::vector<Vec2> const&r,std::vector<Vec2> const&c)->double{
 if(r.empty()||c.size()<2)return std::numeric_limits<double>::max();std::array<std::size_t,4>ix{0,(r.size()-1)/3,2*(r.size()-1)/3,r.size()-1};double sum=0;
 for(auto i:ix){double best=std::numeric_limits<double>::max();for(std::size_t j=0;j+1<c.size();++j){double ax=c[j].x,ay=c[j].y,vx=c[j+1].x-ax,vy=c[j+1].y-ay,vv=vx*vx+vy*vy,t=vv>0?std::clamp(((r[i].x-ax)*vx+(r[i].y-ay)*vy)/vv,0.0,1.0):0.0,dx=r[i].x-(ax+t*vx),dy=r[i].y-(ay+t*vy);best=std::min(best,dx*dx+dy*dy);}sum+=best;}return std::isfinite(sum)?sum:std::numeric_limits<double>::max();
}

struct GlobalRobustNllConfig { double nu{2.4}; };
struct GlobalRobustNllResult { double fitness{std::numeric_limits<double>::max()}; double chi2{0}; bool ok{false}; };

template <typename T>
inline auto evaluateGlobalRobustNll(std::vector<T> const& residual, std::vector<T> const& sigma2,
                                     GlobalRobustNllConfig const& cfg = {}) -> GlobalRobustNllResult {
    GlobalRobustNllResult out;
    auto const n = residual.size();
    if (n == 0 || sigma2.empty()) { return out; }
    bool const broadcast = sigma2.size() == 1;
    if (!broadcast && sigma2.size() != n) { return out; }
    double chi2 = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        double const r = static_cast<double>(residual[i]);
        double const s2 = std::max(static_cast<double>(broadcast ? sigma2[0] : sigma2[i]), 1e-12);
        if (!std::isfinite(r)) { return out; }
        chi2 += (r * r) / s2;
    }
    if (!std::isfinite(chi2) || chi2 < 0.0) { return out; }
    double const nu = std::max(cfg.nu, 1e-6);
    double const val = 0.5 * (nu + static_cast<double>(n)) * std::log1p(chi2 / nu);
    if (!std::isfinite(val)) { return out; }
    out.chi2 = chi2;
    out.fitness = val;
    out.ok = true;
    return out;
}
}
