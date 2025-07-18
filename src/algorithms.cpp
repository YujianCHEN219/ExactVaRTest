#include <Rcpp.h>
#include <unordered_map>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>

using namespace Rcpp;

const double EPS = 1e-15;
static inline double safe_log(double x) { return std::log(x > EPS ? x : EPS); }

/* ---------- helpers ---------- */
static inline long long encode_state(int last,int s00,int s10,int s01,int s11,
                                     long long S2,long long S3,long long S4,long long S5)
{
  return static_cast<long long>(last)
  + S2*static_cast<long long>(s00)
  + S3*static_cast<long long>(s10)
  + S4*static_cast<long long>(s01)
  + S5*static_cast<long long>(s11);
}
struct LLhash { std::size_t operator()(long long x) const { return std::hash<long long>()(x); } };

/* ======================================================================
 *  LR_ind
 * ====================================================================*/
struct IndState {
  long long code = 0;
  int last = 0, s00 = 0, s10 = 0, s01 = 0, s11 = 0;
  double prob = 0.0;
};

// [[Rcpp::export]]
List fb_lrind_fastcpp(int n,
                      double alpha = 0.05,
                      double prune_threshold = 1e-15)
{
  if (n < 1) return List::create(_["LR"] = NumericVector(0),
      _["prob"] = NumericVector(0));
  
  long long S2 = 2LL, S3 = S2*(n+1), S4 = S3*(n+1), S5 = S4*(n+1);
  
  std::vector<IndState> cur;
  cur.reserve(2);
  cur.push_back({encode_state(0,0,0,0,0,S2,S3,S4,S5), 0,0,0,0,0, 1-alpha});
  cur.push_back({encode_state(1,0,0,0,0,S2,S3,S4,S5), 1,0,0,0,0, alpha});
  
  for (int t = 1; t < n; ++t) {
    std::size_t keep = 0;
    for (std::size_t i = 0; i < cur.size(); ++i)
      if (cur[i].prob >= prune_threshold) cur[keep++] = cur[i];
      cur.resize(keep);
      if (cur.empty()) break;
      
      std::unordered_map<long long,IndState,LLhash> nxt;
      nxt.reserve(cur.size()*2);
      
      for (const auto &s : cur) {
        {
          IndState ns = s;
          ns.last ? ++ns.s10 : ++ns.s00;
          ns.last = 0;
          ns.prob *= (1.0-alpha);
          if (ns.prob >= prune_threshold) {
            ns.code = encode_state(ns.last,ns.s00,ns.s10,ns.s01,ns.s11,S2,S3,S4,S5);
            auto &slot = nxt[ns.code];
            if (slot.prob == 0.0) slot = ns; else slot.prob += ns.prob;
          }
        }
        {
          IndState ns = s;
          ns.last ? ++ns.s11 : ++ns.s01;
          ns.last = 1;
          ns.prob *= alpha;
          if (ns.prob >= prune_threshold) {
            ns.code = encode_state(ns.last,ns.s00,ns.s10,ns.s01,ns.s11,S2,S3,S4,S5);
            auto &slot = nxt[ns.code];
            if (slot.prob == 0.0) slot = ns; else slot.prob += ns.prob;
          }
        }
      }
      cur.clear();
      cur.reserve(nxt.size());
      for (auto &kv : nxt) cur.push_back(kv.second);
  }
  
  std::size_t keep = 0;
  for (std::size_t i = 0; i < cur.size(); ++i)
    if (cur[i].prob >= prune_threshold) cur[keep++] = cur[i];
    cur.resize(keep);
    if (cur.empty())
      return List::create(_["LR"] = NumericVector(0),
                          _["prob"] = NumericVector(0));
    
    std::unordered_map<double,double> dist;
    dist.reserve(cur.size());
    for (const auto &s : cur) {
      int T0 = s.s00 + s.s10, T1 = s.s01 + s.s11;
      double pHat = n>1 ? static_cast<double>(T1)/(n-1) : 0.0;
      double num  = T0*safe_log(1-pHat) + T1*safe_log(pHat);
      
      int sum01 = s.s00 + s.s01, sum11 = s.s10 + s.s11;
      double pi01 = sum01 ? static_cast<double>(s.s01)/sum01 : 1.0;
      double pi11 = sum11 ? static_cast<double>(s.s11)/sum11 : 1.0;
      double den  = s.s00*safe_log(1-pi01) + s.s01*safe_log(pi01) +
        s.s10*safe_log(1-pi11) + s.s11*safe_log(pi11);
      
      double LR = -2.0*(num - den);
      if (LR < 0 && LR > -1e-12) LR = 0.0;
      dist[LR] += s.prob;
    }
    
    std::vector<std::pair<double,double>> vec(dist.begin(), dist.end());
    std::sort(vec.begin(), vec.end(),
              [](const auto&a,const auto&b){ return a.first<b.first; });
    
    NumericVector LR_out, P_out;
    for (std::size_t i = 0; i < vec.size();) {
      double curLR = vec[i].first, p = 0.0;
      while (i < vec.size() && vec[i].first == curLR) { p += vec[i].second; ++i; }
      LR_out.push_back(curLR);
      P_out.push_back(p);
    }
    P_out = P_out / std::accumulate(P_out.begin(), P_out.end(), 0.0);
    return List::create(_["LR"] = LR_out, _["prob"] = P_out);
}

/* ======================================================================
 *  LR_cc
 * ====================================================================*/
struct CCState {
  long long code = 0;
  int last = 0, c1 = 0, s00 = 0, s10 = 0, s01 = 0, s11 = 0;
  double prob = 0.0;
};

static inline long long enc6(int last,int c1,int s00,int s10,int s01,int s11,int n)
{
  long long step = static_cast<long long>(n)+10LL;
  long long S2=step, S3=S2*step, S4=S3*step, S5=S4*step, S6=S5*step;
  return last + S2*c1 + S3*s00 + S4*s10 + S5*s01 + S6*s11;
}

// [[Rcpp::export]]
List fb_lrcc_fastcpp(int n,
                     double alpha = 0.05,
                     double prune_threshold = 1e-15)
{
  if (n < 1) return List::create(_["LR"] = NumericVector(0),
      _["prob"] = NumericVector(0));
  
  std::unordered_map<long long,CCState,LLhash> cur;
  cur.reserve(2);
  cur[enc6(0,0,0,0,0,0,n)] = {enc6(0,0,0,0,0,0,n),0,0,0,0,0,0,1-alpha};
  cur[enc6(1,1,0,0,0,0,n)] = {enc6(1,1,0,0,0,0,n),1,1,0,0,0,0,alpha};
  
  for (int t = 1; t < n; ++t) {
    std::unordered_map<long long,CCState,LLhash> nxt;
    nxt.reserve(cur.size()*2);
    
    for (const auto &kv : cur) {
      const CCState &s = kv.second;
      if (s.prob < prune_threshold) continue;
      
      {
        CCState ns = s;
        ns.last ? ++ns.s10 : ++ns.s00;
        ns.last = 0;
        ns.prob *= (1-alpha);
        if (ns.prob >= prune_threshold) {
          ns.code = enc6(ns.last,ns.c1,ns.s00,ns.s10,ns.s01,ns.s11,n);
          auto &slot = nxt[ns.code];
          if (slot.prob == 0.0) slot = ns; else slot.prob += ns.prob;
        }
      }
      {
        CCState ns = s;
        ns.last ? ++ns.s11 : ++ns.s01;
        ++ns.c1;
        ns.last = 1;
        ns.prob *= alpha;
        if (ns.prob >= prune_threshold) {
          ns.code = enc6(ns.last,ns.c1,ns.s00,ns.s10,ns.s01,ns.s11,n);
          auto &slot = nxt[ns.code];
          if (slot.prob == 0.0) slot = ns; else slot.prob += ns.prob;
        }
      }
    }
    cur.swap(nxt);
  }
  
  std::unordered_map<double,double> dist;
  dist.reserve(cur.size());
  for (const auto &kv : cur) {
    const CCState &s = kv.second;
    
    double p_   = std::min(std::max(alpha, EPS), 1.0-EPS);
    double phat = static_cast<double>(s.c1)/n;
    double ph_  = std::min(std::max(phat, EPS), 1.0-EPS);
    double LRuc = -2.0*(s.c1*safe_log(p_) + (n-s.c1)*safe_log(1-p_)
                          - s.c1*safe_log(ph_) - (n-s.c1)*safe_log(1-ph_));
    
    int T0 = s.s00 + s.s10, T1 = s.s01 + s.s11;
    double pHat = n>1 ? static_cast<double>(T1)/(n-1) : 0.0;
    double num  = T0*safe_log(1-pHat) + T1*safe_log(pHat);
    
    int sum01 = s.s00 + s.s01, sum11 = s.s10 + s.s11;
    double pi01 = sum01 ? static_cast<double>(s.s01)/sum01 : 1.0;
    double pi11 = sum11 ? static_cast<double>(s.s11)/sum11 : 1.0;
    double den  = s.s00*safe_log(1-pi01) + s.s01*safe_log(pi01) +
      s.s10*safe_log(1-pi11) + s.s11*safe_log(pi11);
    
    double LRind = -2.0*(num - den);
    double LR = LRuc + LRind;
    if (LR < 0 && LR > -1e-12) LR = 0.0;
    dist[LR] += s.prob;
  }
  
  std::vector<std::pair<double,double>> vec(dist.begin(), dist.end());
  std::sort(vec.begin(), vec.end(),
            [](const auto&a,const auto&b){ return a.first<b.first; });
  
  NumericVector LR_out, P_out;
  for (std::size_t i = 0; i < vec.size();) {
    double curLR = vec[i].first, p = 0.0;
    while (i < vec.size() && vec[i].first == curLR) { p += vec[i].second; ++i; }
    LR_out.push_back(curLR);
    P_out.push_back(p);
  }
  P_out = P_out / std::accumulate(P_out.begin(), P_out.end(), 0.0);
  return List::create(_["LR"] = LR_out, _["prob"] = P_out);
}
