struct fraction{
  int numer,denom;
  
  fraction() : numer(0), denom(1) {}
  fraction(int x) : numer(x), denom(1) {}
  fraction(int x,int y){
    assert(y!=0);
    int G = std::gcd(x,y);
    x /= G,y /= G;
    numer=x,denom=y;
    if(denom<0) numer = -numer, denom = -denom;
  }//x/y
  
  fraction operator+(const fraction &o){
    fraction res;
    res.denom = std::lcm(denom,o.denom);
    res.numer = res.denom/denom*numer + res.denom/o.denom*o.numer;
    int G = std::gcd(res.denom,res.numer);
    res.denom /= G,res.numer /= G;
    return res;
  }
  
  fraction operator+(int x){return *this + fraction(x,1);}
  
  fraction operator-(const fraction &o){
    fraction res(o.numer,-o.denom);
    return *this + res;
  }
  
  fraction operator-(int x){return *this - fraction(x,1);}
  fraction operator-(){return fraction(-numer,denom);}
  
  fraction operator*(const fraction &o){return fraction(numer*o.numer,denom*o.denom);}
  fraction operator*(int x){return fraction(numer*x,denom);}
  fraction operator/(const fraction &o){return fraction(numer*o.denom,denom*o.numer);}
  fraction operator/(int x){return fraction(numer,denom*x);}
  
  fraction operator+=(const fraction &o){*this = *this + o; return *this;}
  fraction operator+=(const int &o){*this = *this + o; return *this;}
  fraction operator-=(const fraction &o){*this = *this - o; return *this;}
  fraction operator-=(const int &o){*this = *this - o; return *this;}
  fraction operator*=(const fraction &o){*this = *this * o; return *this;}
  fraction operator*=(const int &o){*this = *this * o; return *this;}
  fraction operator/=(const fraction &o){*this = *this / o; return *this;}
  fraction operator/=(const int &o){*this = *this / o; return *this;}
  
  fraction operator%(int x){
    numer %= (x*denom);
    if(numer<0) numer += x*denom;
    return fraction(numer,denom);
  }
  
  
  bool operator==(const fraction &o) const{return numer*o.denom == o.numer*denom;}
  bool operator!=(const fraction &o) const{return !(*this == o);}
  bool operator==(int x) const{return numer == x*denom;}
  bool operator!=(int x) const{return !(*this==x);}
  bool operator<(const fraction &o) const{return numer*o.denom < o.numer*denom;}
  bool operator>(const fraction &o) const{return o < *this;}
  bool operator<=(const fraction &o) const{return !(o < *this);}
  bool operator>=(const fraction &o) const{return !(o > *this);}
  
  friend std::ostream& operator<<(std::ostream& os,const fraction &x){
    os<<x.numer<<'/'<<x.denom;
    return os;
  }
};
