#include <bits/stdc++.h>
#include <omp.h>

using namespace std;

typedef tuple<int, int> point;
typedef tuple<point, point> line;
typedef tuple<point, point, point> triangle;

int nextInt() {
  int x;
  scanf("%d", &x);
  return x;
}

int N;
point s, e;

vector<triangle> obst;
vector<line> L;
vector<point> P;
map< point, vector<point> > G;

inline int minX(point a, point b) { return min(get<0>(a), get<0>(b)); }
inline int minY(point a, point b) { return min(get<1>(a), get<1>(b)); }
inline int maxX(point a, point b) { return max(get<0>(a), get<0>(b)); }
inline int maxY(point a, point b) { return max(get<1>(a), get<1>(b)); }

double distance(point a, point b) {
  return sqrt(pow(1. * (get<0>(a) - get<0>(b)), 2) +
              pow(1. * (get<1>(a) - get<1>(b)), 2));
}

bool onSegment(point p, point q, point r) {
  if (get<0>(q) <= maxX(p, r) && get<0>(q) >= minX(p, r) &&
      get<1>(q) <= maxY(p, r) && get<1>(q) >= minY(p, r))
    return true;
  return false;
}

int orientation(point p, point q, point r) {
  int val = (get<1>(q) - get<1>(p)) * (get<0>(r) - get<0>(q)) -
            (get<0>(q) - get<0>(p)) * (get<1>(r) - get<1>(q));
  if (val == 0) return 0;
  return (val > 0) ? 1 : 2;
}

bool intersect(point p1, point q1, point p2, point q2) {
  int o1 = orientation(p1, q1, p2);
  int o2 = orientation(p1, q1, q2);
  int o3 = orientation(p2, q2, p1);
  int o4 = orientation(p2, q2, q1);
  if (o1 != o2 && o3 != o4) return true;
  if (o1 == 0 && onSegment(p1, p2, q1)) return true;
  if (o2 == 0 && onSegment(p1, q2, q1)) return true;
  if (o3 == 0 && onSegment(p2, p1, q2)) return true;
  if (o4 == 0 && onSegment(p2, q1, q2)) return true;
  return false;
}

double sign (point p1, point p2, point p3)
{
    return (get<0>(p1) - get<0>(p3)) * (get<1>(p2) - get<1>(p3)) - (get<0>(p2) - get<0>(p3)) * (get<1>(p1) - get<1>(p3));
}

bool PointInTriangle (point pt, point v1, point v2, point v3)
{
    bool b1 = sign(pt, v1, v2) < 0.0;
    bool b2 = sign(pt, v2, v3) < 0.0;
    bool b3 = sign(pt, v3, v1) < 0.0;
    return ((b1 == b2) && (b2 == b3));
}

int main() {
  int a0 = nextInt();
  int a1 = nextInt();
  int a2 = nextInt();
  int a3 = nextInt();

  s = make_tuple(a0, a1);
  e = make_tuple(a2, a3);

  N = nextInt();
  obst.reserve(N);

  for (int i = 0; i < N; i++) {
    int x1 = nextInt();
    int x2 = nextInt();
    int x3 = nextInt();
    int x4 = nextInt();
    int x5 = nextInt();
    int x6 = nextInt();

    obst.push_back(
        make_tuple(make_tuple(x1, x2), make_tuple(x3, x4), make_tuple(x5, x6)));
  }

  //L.reserve(3 * N);
  //P.reserve(N * 12 + 2);

  P.push_back(s);
  P.push_back(e);

  for (auto o : obst) {
    auto p0 = get<0>(o);
    auto p1 = get<1>(o);
    auto p2 = get<2>(o);

    // printf("%d %d %d %d %d %d\n", get<0>(p0), get<1>(p0), get<0>(p1),
    // get<1>(p1), get<0>(p2), get<1>(p2));
    
    auto l01 = get<0>(p0) < get<0>(p1) ? make_tuple(p0, p1) : make_tuple(p1, p0);
    auto l02 = get<0>(p0) < get<0>(p2) ? make_tuple(p0, p2) : make_tuple(p2, p0);
    auto l12 = get<0>(p1) < get<0>(p2) ? make_tuple(p1, p2) : make_tuple(p2, p1);

    auto p00 = make_tuple(get<0>(p0) + 1, get<1>(p0));
    auto p01 = make_tuple(get<0>(p0), get<1>(p0) + 1);
    auto p02 = make_tuple(get<0>(p0) - 1, get<1>(p0));
    auto p03 = make_tuple(get<0>(p0), get<1>(p0) - 1);

    auto p10 = make_tuple(get<0>(p1) + 1, get<1>(p1));
    auto p11 = make_tuple(get<0>(p1), get<1>(p1) + 1);
    auto p12 = make_tuple(get<0>(p1) - 1, get<1>(p1));
    auto p13 = make_tuple(get<0>(p1), get<1>(p1) - 1);

    auto p20 = make_tuple(get<0>(p2) + 1, get<1>(p2));
    auto p21 = make_tuple(get<0>(p2), get<1>(p2) + 1);
    auto p22 = make_tuple(get<0>(p2) - 1, get<1>(p2));
    auto p23 = make_tuple(get<0>(p2), get<1>(p2) - 1);

    L.push_back(l01);
    L.push_back(l02);
    L.push_back(l12);

    if( !PointInTriangle(p00, p0, p1, p2) ) P.push_back(p00);
    if( !PointInTriangle(p01, p0, p1, p2) ) P.push_back(p01);
    if( !PointInTriangle(p02, p0, p1, p2) ) P.push_back(p02);
    if( !PointInTriangle(p03, p0, p1, p2) ) P.push_back(p03);

    if( !PointInTriangle(p10, p0, p1, p2) ) P.push_back(p10);
    if( !PointInTriangle(p11, p0, p1, p2) ) P.push_back(p11);
    if( !PointInTriangle(p12, p0, p1, p2) ) P.push_back(p12);
    if( !PointInTriangle(p13, p0, p1, p2) ) P.push_back(p13);

    if( !PointInTriangle(p20, p0, p1, p2) ) P.push_back(p20);
    if( !PointInTriangle(p21, p0, p1, p2) ) P.push_back(p21);
    if( !PointInTriangle(p22, p0, p1, p2) ) P.push_back(p22);
    if( !PointInTriangle(p23, p0, p1, p2) ) P.push_back(p23);

    /*
    auto dp00 = make_tuple(get<0>(p0) + 1, get<1>(p0) + 1);
    auto dp01 = make_tuple(get<0>(p0) + 1, get<1>(p0) - 1);
    auto dp02 = make_tuple(get<0>(p0) - 1, get<1>(p0) + 1);
    auto dp03 = make_tuple(get<0>(p0) - 1, get<1>(p0) - 1);

    auto dp10 = make_tuple(get<0>(p1) + 1, get<1>(p1) + 1);
    auto dp11 = make_tuple(get<0>(p1) + 1, get<1>(p1) - 1);
    auto dp12 = make_tuple(get<0>(p1) - 1, get<1>(p1) + 1);
    auto dp13 = make_tuple(get<0>(p1) - 1, get<1>(p1) - 1);

    auto dp20 = make_tuple(get<0>(p2) + 1, get<1>(p2) + 1);
    auto dp21 = make_tuple(get<0>(p2) + 1, get<1>(p2) - 1);
    auto dp22 = make_tuple(get<0>(p2) - 1, get<1>(p2) + 1);
    auto dp23 = make_tuple(get<0>(p2) - 1, get<1>(p2) - 1);

    P.push_back(dp00);
    P.push_back(dp01);
    P.push_back(dp02);
    P.push_back(dp03);

    P.push_back(dp10);
    P.push_back(dp11);
    P.push_back(dp12);
    P.push_back(dp13);

    P.push_back(dp20);
    P.push_back(dp21);
    P.push_back(dp22);
    P.push_back(dp23);
    */
  }

  printf("P prima %d\n", P.size());
  vector<point> tmp;
  #pragma omp parallel for schedule(guided)
  for(int i=0; i<P.size(); ++i)
  {
    auto p = P[i];
    bool out = true;
    for(auto o : obst)
    {
      if( PointInTriangle(p, get<0>(o), get<1>(o), get<2>(o)) ) 
      {
        out = false;
        break;
      }
    }
    #pragma omp critical 
    {
      if( out ) tmp.push_back(p);
    }
  }
  P = tmp;
  printf("P dopo %d\n", P.size());
  printf("L dopo %d\n", L.size());
  fflush(stdout);

  auto linecp = [](line a, line b){ return get<0>(get<0>(a)) < get<0>(get<0>(b)); };

  sort(P.begin(), P.end());
  sort(L.begin(), L.end(), linecp);
  printf("SORT\n"); 
  fflush(stdout);

/*
  #pragma omp parallel for schedule(guided)
  for(int i=0; i<P.size(); i++)
  {
    for(int j=0; j<P.size(); j++)
    {
      if( i == j ) continue;
      bool out = true;

      for (auto l : L) {
        if( get<0>(get<0>(l)) > get<0>(P[i]) ) break;
        if (intersect(P[i], P[j], get<0>(l), get<1>(l))) {
          out = false;
          break;
        }
      }
      if( out ){
        G[P[i]].push_back(P[j]);
      }
    }
  }

  int cc = 0;

  for(auto g : G)
    cc += g.second.size();

  printf("|G| = %d\n", cc);
  fflush(stdout);
*/

//   return 0;

//  for( auto l : L )
//  {
//    printf("[%d][%d] - [%d][%d]\n", get<0>(get<0>(l)),get<1>(get<0>(l)),get<0>(get<1>(l)),get<1>(get<1>(l))); 
//  }

  set<point> S(P.begin(), P.end());
  map<point, double> best;
  map<point, point> prec;
  priority_queue<pair<double, point>> Q;

  Q.push(make_pair(0., s));
  best[s] = 0.;
  best[e] = INT_MAX;

  int count = 0;
  while (!Q.empty()) {
    auto top = Q.top();
    Q.pop();

    if (best[top.second] != -top.first) continue;
    if (S.find(top.second) == S.end()) continue;
    S.erase(top.second);

    printf("[%d][%d][%.6f][%d]\n", get<0>(top.second), get<1>(top.second), top.first, count++);
    fflush(stdout);
    // printf("[%d, %d] [%.6f]\n", get<0>(top.second), get<1>(top.second),
    // best[top.second]);
    if (top.second == e) break;  // END

    #pragma omp parallel for schedule(guided)
    for (int i = 0; i < P.size(); ++i) {

      auto p = P[i];

      if (S.find(p) == S.end()) continue;
      if( abs(get<0>(p) - get<0>(top.second)) > 1000 ) continue;
      if( abs(get<1>(p) - get<1>(top.second)) > 1000 ) continue;

      bool inter = false;
      auto pp = make_tuple(top.second, p);
      if( get<0>(top.second) > get<0>(p) )
        pp = make_tuple(p, top.second);

      auto begin = lower_bound(L.begin(), L.end(), pp, linecp);
      auto end = L.end(); // upper_bound(L.begin(), L.end(), pp, linecp)+1;

      for(begin; begin != end; begin++)
      {
        auto l = *begin;
//        if( get<0>(get<1>(l)) > get<0>(get<0>(pp)) ) continue;
        if( get<0>(get<0>(l)) > get<0>(get<1>(pp)) ) break;
        if( intersect(top.second, p, get<0>(l), get<1>(l)) )
        {
          inter = true;
          break;
        }
      }

      if( inter ) continue;

      begin = lower_bound(L.begin(), L.end(), pp, linecp);
      end = L.begin(); // upper_bound(L.begin(), L.end(), pp, linecp)+1;

      for(begin; begin != end; begin--)
      {
        auto l = *begin;
//        if( get<0>(get<1>(l)) > get<0>(get<0>(pp)) ) continue;
        if( intersect(top.second, p, get<0>(l), get<1>(l)) )
        {
          inter = true;
          break;
        }
      }

      if( inter ) continue;

      double dist = best[top.second] + distance(top.second, p);
      if (best.find(p) == best.end() || dist < best[p]) {
        #pragma omp critical
        {
          best[p] = dist;
          prec[p] = top.second;
          Q.push(make_pair(-dist, p));
        }
      }
    }
    

    /*
    #pragma omp parallel for schedule(guided)
    for (int i = 0; i < G[top.second].size(); ++i) {

      auto p = G[top.second][i];

      if (S.find(p) == S.end()) continue;

      double dist = best[top.second] + distance(top.second, p);
      if (best.find(p) == best.end() || dist < best[p]) {
        #pragma omp critical
        {
          best[p] = dist;
          prec[p] = top.second;
          Q.push(make_pair(-dist, p));
        }
      }
    }
    */
  }

  if (prec.find(e) == prec.end()) {
    printf("IMPOSSIBLE");
  } else {
    vector<point> output;
    point last = e;
    do {
      output.push_back(last);
      last = prec[last];
    } while (last != s);
    output.push_back(s);
    reverse(output.begin(), output.end());

    printf("%d\n", output.size());
    for (auto o : output) {
      printf("%d %d\n", get<0>(o), get<1>(o));
    }
  }
  fflush(stdout);
  return 0;
}
