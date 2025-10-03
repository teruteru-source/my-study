//コードはC++で記述

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cassert>

using D = double;
const D EPS = 1e-10;//丸め誤差を考える

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
        T ma = EPS;
        for (int row = rank; row < m; ++row) {
            if (abs(A[row][col]) > ma) {
                ma = abs(A[row][col]);
                pivot = row;
            }
        }
        if (pivot == -1) continue;

        swap(A[pivot], A[rank]);

        auto fac = A[rank][col];
        for (int col2 = 0; col2 < n; ++col2) A[rank][col2] /= fac;

        for (int row = 0; row < m; ++row) {
            if (row != rank && abs(A[row][col]) > EPS) {
                auto fac = A[row][col];
                for (int col2 = 0; col2 < n; ++col2) {
                    A[row][col2] -= A[rank][col2] * fac;
                }
            }
        }
        ++rank;
    }
    return rank;
}

template<class T> std::pair<std::vector<int>,std::vector<std::vector<double>>> linear_equation(int p,Matrix<T> A, std::vector<T> b) {
  //連立1次方程式を解いて基底および関係式を返す
    int m = A.size(), n = A[0].size();
    Matrix<T> M(m, n + 1);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) M[i][j] = A[i][j];
        M[i][n] = b[i];
    }
    int rank = GaussJordan(M, true);
    std::cout<<"The dimension of H("<<p<<") is: "<<n-rank<<"\n";
    std::vector<int> basis;
    int j=0;
    for(int i=0;i<m;i++){
      if(j>=n) break;
      while(j<n && M[i][j] < EPS){
        basis.push_back(j);
        j++;
      }
      j++;
    }
    //for (int row = rank; row < m; ++row) if (abs(M[row][n]) > EPS) return {-11111};
    std::cout<<"The bases are: ";
    for(int i=0;i<basis.size();i++){
      std::cout<<"(";
      if(basis[i]<n-1) std::cout<<basis[i];
      else std::cout<<"inf";
      std::cout<<")";
      if(i<int(basis.size())-1) std::cout<<",";
      else std::cout<<"\n";
    }
    
    std::vector<std::vector<double>> ans(m,std::vector<double>(basis.size(),0));
    j=0;
    for(int i=0;i<M.size();i++){
      while(j<n && M[i][j] < EPS) j++;
      if(j>=n) break;
      
      std::cout<<"("<<j<<") = ";
      if(std::count(M[i].begin(),M[i].end(),0) == n) std::cout<<"0\n";
      else{
        for(int k=0;k<basis.size();k++) if(abs(M[i][basis[k]]) > EPS){
          if(j!=basis[k]) ans[i][k] = -M[i][basis[k]];
          else ans[i][k] = M[i][basis[k]];
          std::cout<<-M[i][basis[k]]<<" * (";
          if(basis[k]<n-1) std::cout<<basis[k];
          else std::cout<<"inf";
          std::cout<<"),";
        }
      }
    }
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

std::vector<Matrix<int>> Heilbronn_Matrix(int p){
  //素数pでのHeilbronn Matrixの集合を返す
  //ここらへんがバグっている？
  assert(is_Prime(p));
  std::vector<Matrix<int>> result;
  Matrix<int> init(2,2);
  init.val[0][0]=1,init.val[0][1]=0,init.val[1][0]=0,init.val[1][1]=p;
  result.push_back(init);
	for(int r=-(p-1)/2;r<=p/2;r++)
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

std::vector<double> Calc_Hecke_operator(int p,std::vector<int> &basis,std::vector<std::vector<double>> &state,int s){
  //素数pでのHecke作用素の値を得る,sは代入するM-symbols
  std::vector<Matrix<int>> Heilbronn = Heilbronn_Matrix(p);
  std::vector<double> result(basis.size(),0);
  for(auto c : Heilbronn){
    std::cerr<<c.val[0][0]<<" "<<c.val[0][1]<<" "<<c.val[1][0]<<" "<<c.val[1][1]<<"\n";
    int m1 = (s*c.val[0][0] + c.val[1][0]) % p;
    int m2 = (s*c.val[0][1] + c.val[1][1]) % p;
    int num = convert_in_prime_forms(M_symbols{m1,m2},p);
    for(int i=0;i<basis.size();i++) result[i] += state[num][i];
  }
  return result;
}

int main(){
  int p;
  std::cout<<"Input the prime p for which you want to know L(f,1)/Omega(f) of the newform f of level p:\n";
  std::cin>>p;
  if(!is_Prime(p)){
    std::cout<<p<<" is not a prime :(\n";
    return 1;
  }
  Matrix<double> mat(2*p,p+1);
  for(int i=0;i<p;i++){//二項関係式
    mat.val[i][i] = 1;
    mat.val[i][convert_in_prime_forms(M_symbols{-1,i},p)] = 1;
  }
  
  for(int i=0;i<p;i++){//三項関係式
    mat.val[i+p][i] = 1;
    mat.val[i+p][convert_in_prime_forms(M_symbols{i+1,-i},p)] = 1;
    mat.val[i+p][convert_in_prime_forms(M_symbols{1,-i-1},p)] = 1;
  }
  std::vector<double> zeros(2*p,0);
  auto [basis,state] = linear_equation(p,mat,zeros);
  auto result = Calc_Hecke_operator(2,basis,state,basis[0]);
  for(int i=0;i<basis.size();i++) std::cout<<result[i]<<" ";
  std::cout<<"L(f,1)/Omega(f) = 1/"<<(3-result[0])<<"\n";
  return 0;
}
