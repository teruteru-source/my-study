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
