//コードはC++で記述

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <numeric>
#include <iomanip>

std::vector<int> Parametrization_Degrees = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,-1,-1,1,1,-1,1,-1,1,1,1,-1,-1,1,-1,2,1,-1,-1,2};

struct fraction{
  int numerator,denominator;
  
  fraction() : numerator(0), denominator(1) {}
  fraction(int x) : numerator(x), denominator(1) {}
  fraction(int x,int y){
    assert(y!=0);
    int G = std::gcd(x,y);
    x /= G,y /= G;
    numerator=x,denominator=y;
    if(denominator<0) numerator = -numerator, denominator = -denominator;
  }//x/y
  
  fraction operator+(const fraction &o){
    fraction res;
    res.denominator = std::lcm(denominator,o.denominator);
    res.numerator = res.denominator/denominator*numerator + res.denominator/o.denominator*o.numerator;
    int G = std::gcd(res.denominator,res.numerator);
    res.denominator /= G,res.numerator /= G;
    return res;
  }
  fraction operator+(int x){return *this + fraction(x,1);}
  
  fraction operator-(const fraction &o){
    fraction res(o.numerator,-o.denominator);
    return *this + res;
  }
  fraction operator-(int x){return *this - fraction(x,1);}
  fraction operator-(){return fraction(-numerator,denominator);}
  
  fraction operator*(const fraction &o){return fraction(numerator*o.numerator,denominator*o.denominator);}
  fraction operator*(int x){return fraction(numerator*x,denominator);}
  fraction operator/(const fraction &o){return fraction(numerator*o.denominator,denominator*o.numerator);}
  fraction operator/(int x){return fraction(numerator,denominator*x);}
  
  fraction operator+=(const fraction &o){*this = *this + o; return *this;}
  fraction operator+=(const int &o){*this = *this + o; return *this;}
  fraction operator-=(const fraction &o){*this = *this - o; return *this;}
  fraction operator-=(const int &o){*this = *this - o; return *this;}
  fraction operator*=(const fraction &o){*this = *this * o; return *this;}
  fraction operator*=(const int &o){*this = *this * o; return *this;}
  fraction operator/=(const fraction &o){*this = *this / o; return *this;}
  fraction operator/=(const int &o){*this = *this / o; return *this;}
  
  fraction operator%(int x){
    numerator %= (x*denominator);
    if(numerator<0) numerator += x*denominator;
    return fraction(numerator,denominator);
  }
  
  bool operator==(const fraction &o) const{return numerator*o.denominator == o.numerator*denominator;}
  bool operator!=(const fraction &o) const{return !(*this == o);}
  bool operator==(int x) const{return numerator == x*denominator;}
  bool operator!=(int x) const{return !(*this==x);}
  bool operator<(const fraction &o) const{return numerator*o.denominator < o.numerator*denominator;}
  bool operator>(const fraction &o) const{return o < *this;}
  bool operator<=(const fraction &o) const{return !(o < *this);}
  bool operator>=(const fraction &o) const{return !(o > *this);}
  
  friend std::ostream& operator<<(std::ostream& os,const fraction &x){
    os<<x.numerator<<'/'<<x.denominator;
    return os;
  }
  
  fraction inv(){
    assert(numerator != 0);
    return fraction(denominator,numerator);
  }
  
  std::vector<int> continued_fraction(){//連分数展開
    std::vector<int> res;
    int p=numerator,q=denominator;
    while(q != 0){
      int a = p/q, r = p%q;
      res.push_back(a);
      p=q,q=r;
    }
    return res;
  }
  
  std::vector<fraction> Construct_Fractional_Convergence_Sequence(){//分数収束列の構成
    std::vector<int> C = this->continued_fraction();
    //for(int i=0;i<C.size();i++) cerr<<C[i]<<"\n";
    std::vector<fraction> res;
    int pre_p = 1, pre_q = 0, p = C[0], q = 1;
    res.push_back(fraction(p,q));
    for(int i=1;i<C.size();i++){
      int next_p = C[i]*p + pre_p;
      int next_q = C[i]*q + pre_q;
      pre_p=p,pre_q=q,p=next_p,q=next_q;
      res.push_back(fraction(p,q));
    }
    return res;
  }
  
  double val(){
    return (double)numerator / (double)denominator;
  }
};

fraction abs(fraction &x){return std::max(-x,x);}

template<class T> struct Matrix {//行列の構造体
    std::vector<std::vector<T> > val;
    Matrix(int n, int m, T x = 0) : val(n, std::vector<T>(m, x)) {}
    void init(int n, int m, T x = 0) {val.assign(n, std::vector<T>(m, x));}
    size_t size() const {return val.size();}
    inline std::vector<T>& operator [] (int i) {return val[i];}
};

template<class T> int GaussJordan(Matrix<T> &A, bool is_extended = false) {
  //掃き出し法
  //連立1次方程式を解くときは，is_extended = true
    int m = A.size(), n = A[0].size();
    int rank = 0;
    for (int col = 0; col < n; ++col) {
        if (is_extended && col == n-1) break;

        int pivot = -1;
        T ma = T(0);
        for (int row = rank; row < m; ++row) {
            if (abs(A[row][col]) > ma) {
                ma = A[row][col];
                pivot = row;
            }
        }
        if (pivot == -1) continue;

        swap(A[pivot], A[rank]);

        T fac = A[rank][col];
        for (int col2 = 0; col2 < n; ++col2) A[rank][col2] /= fac;

        for (int row = 0; row < m; ++row) {
            if (row != rank && A[row][col] != 0) {
                fac = A[row][col];
                for (int col2 = 0; col2 < n; ++col2) {
                    A[row][col2] -= A[rank][col2] * fac;
                }
            }
        }
        ++rank;
      /*for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
          std::cerr<<A[i][j]<<" ";
        }
      std::cerr<<"\n";
    }*/
  }
  return rank;
}

template<class T> std::pair<std::vector<int>,std::vector<std::vector<T>>> linear_equation(Matrix<T> A, std::vector<T> b) {
  //連立1次方程式を解いて基底および関係式を返す
    int m = A.size(), n = A[0].size();
    Matrix<T> M(m, n + 1);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) M[i][j] = A[i][j];
        M[i][n] = b[i];
    }
    int rank = GaussJordan(M, true);
    /*for(int i=0;i<m;i++){
      for(int j=0;j<=n;j++){
        std::cerr<<M[i][j]<<" ";
      }
      std::cerr<<"\n";
    }*/
    //std::cout<<"The dimension of H("<<p<<") is: "<<n-rank<<"\n";
    std::vector<int> basis;
    int j=0;
    for(int i=0;i<m;i++){
      if(j>=n) break;
      while(j<n && M[i][j] == 0){
        basis.push_back(j);
        j++;
      }
      j++;
    }
    /*for(auto c : basis) std::cerr<<c<<" ";
    std::cerr<<"\n";*/
    std::vector<std::vector<T>> ans(n+1,std::vector<T>(basis.size(),0));
    for(int i=0;i<m;i++){
      if(std::count(M[i].begin(),M[i].end(),T(0)) == n+1) break;
      else{
        int B = -1;
        for(int j=0;j<=n;j++) if(M[i][j] != T(0) && std::count(basis.begin(),basis.end(),j) == 0) B = j;
        for(int k=0;k<int(basis.size());k++){
          ans[B][k] = -M[i][basis[k]] / M[i][B];
        }
      }
    }
    for(int i=0;i<basis.size();i++) ans[basis[i]][i] = T(1);
    /*for(int i=0;i<n+1;i++){
      for(int j=0;j<basis.size();j++){
        std::cerr<<ans[i][j]<<" ";
      }
      std::cerr<<"\n";
    }*/
    return {basis,ans};
}

struct M_symbols{int c,d;};//M-symbolsの構造体

int mod_pow(int x,int y,int p){//x^y mod p
  x = ((x%p)+p)%p;
  int res = 1;
  while(y>0){
    if(y&1) res *= x;
    x *= x;
    res %= p;
    x %= p;
    y >>= 1;
  }
  return res;
}

int convert_in_prime_forms(M_symbols M,int p){//M-symbolsを(a)の形に変換する
  if(M.d%p==0) return p;//(∞)に対応
  else return (M.c * mod_pow(M.d,p-2,p)+p)%p;
}

M_symbols convert_in_normal_forms(int x,int p){//(a)の形のM-Symbolsを普通の形にする
  if(x==p) return M_symbols{1,0};
  else return M_symbols{x,1};
}

bool is_Prime(int p){//pが素数かどうか判定
  for(int i=2;i*i<=p;i++) if(p % i == 0) return false;
  return true;
}

struct Modular_Symbols{fraction num1,num2;};

M_symbols Modular_Symbol_to_M_Symbols(Modular_Symbols M){
  return M_symbols{(M.num1.denominator * M.num2.numerator - M.num1.numerator * M.num2.denominator)*M.num2.denominator,M.num1.denominator};
}

std::vector<int> Modular_Symbols_to_M_Symbols(Modular_Symbols M,int p){
  assert(is_Prime(p));
  std::vector<fraction> res1 = M.num1.Construct_Fractional_Convergence_Sequence();
  std::vector<fraction> res2 = M.num2.Construct_Fractional_Convergence_Sequence();
  //for(auto c : res2) std::cerr<<c<<"\n";
  int l1 = res1.size(), l2 = res2.size();
  std::vector<int> cnt(p+1,0);
  //for(int i=0;i<l2-1;i++) cerr<<Modular_Symbol_to_M_Symbols(Modular_Symbols{res2[i],res2[i+1]},p)<<"\n";
  for(int i=0;i<l2-1;i++){
    M_symbols ms = Modular_Symbol_to_M_Symbols(Modular_Symbols{res2[i],res2[i+1]});
    cnt[convert_in_prime_forms(ms,p)]++;
  }
  for(int i=0;i<l1-1;i++){
    M_symbols ms = Modular_Symbol_to_M_Symbols(Modular_Symbols{res1[i],res1[i+1]});
    cnt[convert_in_prime_forms(ms,p)]--;
  }
  return cnt;
}

int Legendre(int a,int p){//Legendre記号をEuler規準で計算
  assert(is_Prime(p) && p>2);
  if(a%p == 0) return 0;
  else{
    int res = mod_pow(a,(p-1)/2,p);
    if(res==p-1) res -= p;
    return res;
  }
}

std::vector<int> gamma_l(int p,int L){
  std::vector<int> cnt(p+1,0);
  for(int l=1;l<L;l++){
    std::vector<int> pcnt = Modular_Symbols_to_M_Symbols(Modular_Symbols{fraction(0),fraction(l,L)},p);
    for(int i=0;i<=p;i++) cnt[i] += Legendre(l,p) * pcnt[i];
  }
  return cnt;
}

std::vector<Matrix<int>> Heilbronn_Matrix(int p){
  //素数pでのHeilbronn Matrixの集合を返す
  assert(is_Prime(p));
  std::vector<Matrix<int>> result;
  Matrix<int> init(2,2);
  init.val[0][0]=1,init.val[0][1]=0,init.val[1][0]=0,init.val[1][1]=p;
  result.push_back(init);
	for(int r=-p/2;r<(p+1)/2;r++)
	{
		int x1=p,x2=-r,y1=0,y2=1,a=-p,b=r;
		Matrix<int> now(2,2);
		now.val[0][0] = x1;
		now.val[0][1] = x2;
		now.val[1][0] = y1;
		now.val[1][1] = y2;
		result.push_back(now);
		while(b != 0){
			//if(a<0 && b<0) {a=-a;b=-b;}
			int q = std::round(double(a)/double(b));//最も近い整数
			int c = a-b*q;
			a=-b;
			b=c;
			int x3 = q*x2-x1;
			x1=x2,x2=x3;
			int y3 = q*y2-y1;
			y1=y2,y2=y3;
			now.val[0][0] = x1;
		  now.val[0][1] = x2;
		  now.val[1][0] = y1;
		  now.val[1][1] = y2;
		  result.push_back(now);
		}
	}
	return result;
}

std::vector<fraction> Calc_Hecke_operator(int p,int q,std::vector<int> &basis,std::vector<std::vector<fraction>> &state,int s){
  //素数pでのHecke作用素の値を得る,sは代入するM-symbols
  std::vector<Matrix<int>> Heilbronn = Heilbronn_Matrix(p);
  std::vector<fraction> result(basis.size(),0);
  for(auto c : Heilbronn){
    //std::cerr<<c.val[0][0]<<" "<<c.val[0][1]<<" "<<c.val[1][0]<<" "<<c.val[1][1]<<"\n";
    int m1 = (s*c.val[0][0] + c.val[1][0]) % q;
    int m2 = (s*c.val[0][1] + c.val[1][1]) % q;
    int num = convert_in_prime_forms(M_symbols{m1,m2},q);
    //std::cerr<<num<<"\n";
    for(int i=0;i<int(basis.size());i++) result[i] += state[num][i];
  }
  return result;
}

bool is_rectangular(std::vector<int> &basis,std::vector<std::vector<fraction>> &state,std::vector<fraction> &v_plus, std::vector<fraction> &v_minus,int p){
  //単位格子がrectangularかどうか返す
  assert(state[0].size() == 3);
  v_plus[1] = state[p-basis[0]][1]+1;
  v_plus[0] = -state[p-basis[1]][1];
  v_minus[1] = state[p-basis[0]][1]-1;
  v_minus[0] = -state[p-basis[1]][1];
  return !(((v_plus[0]-v_minus[0])%2) == 0 && ((v_plus[1]-v_minus[1])%2) == 0);
}

int main(){
  std::cout<<
  "Please Input the prime numbers p and l for which you wish to know the value, L(sym^2(f),2)/(L(f,1)*L(f TENSOR chi_l,1)). :\n";
  int p,l;
  std::cin>>p>>l;
  if(!is_Prime(p)){
    std::cout<<"p is not a prime :(\n";
    return 1;
  }
  if(!is_Prime(l) || l%4!=3 || l%p == 0){
    std::cout<<"l is not valid :(\n";
    return 1;
  }
  Matrix<fraction> mat(2*p,p+1);
  for(int i=0;i<p;i++){//二項関係式
    //std::cout<<i<<" "<<convert_in_prime_forms(M_symbols{-1,i},p)<<"\n";
    mat.val[i][i] = fraction(1);
    mat.val[i][convert_in_prime_forms(M_symbols{-1,i},p)] = fraction(1);
  }
  
  for(int i=0;i<p;i++){//三項関係式
    //std::cout<<i<<" "<<convert_in_prime_forms(M_symbols{i+1,-i},p)<<" "<<convert_in_prime_forms(M_symbols{1,-i-1},p)<<"\n";
    mat.val[i+p][i] = fraction(1);
    mat.val[i+p][convert_in_prime_forms(M_symbols{i+1,-i},p)] = fraction(1);
    mat.val[i+p][convert_in_prime_forms(M_symbols{1,-i-1},p)] = fraction(1);
  }
  
  std::vector<fraction> zeros(2*p,0);
  auto [basis,state] = linear_equation(mat,zeros);
  std::vector<fraction> v_plus(2),v_minus(2);
  bool rectangular = is_rectangular(basis,state,v_plus,v_minus,p);
  std::cout<<v_minus[0]<<" "<<v_minus[1]<<"\n";
  std::cout<<v_plus[0]<<" "<<v_plus[1]<<"\n";
  if(rectangular) std::cout<<"The period lattice is rectangular\n";
  else std::cout<<"The period lattice is non-rectangular\n";
  auto res = Calc_Hecke_operator(2,p,basis,state,basis[0]);
  assert(res[0].denominator == 1);
  fraction L_f1_over_Omega_f = fraction(1,(3-res[0].numerator));
  //==1/(1+p-a_p(f)),p=2より
  std::cout<<"L(f,1)/Omega(f) = "<<L_f1_over_Omega_f<<"\n";
  
  std::vector<int> cnt = gamma_l(p,l);
  for(int i=0;i<=p;i++) std::cout<<cnt[i]<<" ";
  std::cout<<"\n";
  std::vector<fraction> bases(basis.size(),0);
  for(int i=0;i<=p;i++) for(int j=0;j<basis.size();j++) bases[j] += (state[i][j] * cnt[i]);
  for(int i=0;i<basis.size();i++) std::cout<<bases[i]<<" ";
  std::cout<<"\n";
  fraction m_minus = 0;
  for(int i=0;i<ssize(basis)-1;i++) m_minus += v_minus[i] * bases[i];
  std::cout<<"m_minus = "<<m_minus<<"\n";
  if(m_minus == 0){
    std::cout<<"m_minus = 0";
    return 0;
  }
  fraction result = L_f1_over_Omega_f.inv()/2 * fraction(1,p) * Parametrization_Degrees[p] * m_minus.inv();
  if(rectangular) result *= 2;
  else result *= 4;
  
  std::cout<<"L(sym^2(f),2)/(L(f,1)*L(f TENSOR chi_l,1)) = "<<result<<"*pi*sqrt("<<l<<")\n";
  double value = result.val() * std::numbers::pi * sqrt(l);
  std::cout<<std::fixed<<std::setprecision(10)<<"= "<<value<<"\n";
  return 0;
}
