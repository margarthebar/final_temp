#include <stdio.h>
#include <stdlib.h>

int order_points(int d01, int d02, int d12){
  int ans;
  if(d01>0){
    if(d02<0){
      //+--:201
      ans = 201;
    }else{
      if(d12>0){
	//+++:012
	ans = 12;
      }else{
	//++-:021
	//++0:021
	//+0-:021
	ans = 21;
      }
    }
  }else{
    if(d02>0){
      //-++:102
      //0++:102
      ans = 102;
    }else{
      if(d12<0){
	//---:210
	//0--:210
        ans = 210;
      }else{
	//--+:120
	//-0+:120
	//--0:120
	ans = 120;
      }
    }
  }
  return ans;
}

int main(){
  int x0,y0,x1,y1,x2,y2,d01,d02,d12,T,M,B,ans;
  x0 = 4;
  y0 = 2;
  x1 = 3;
  y1 = 4;
  x2 = 1;
  y2 = 1;
  d01 = y0 - y1;
  d02 = y0 - y2;
  d12 = y1 - y2;
  ans = order_points(d01, d02, d12);
  T = ans/100;
  M = (ans%100)/10;
  B = (ans%10);
  printf("T:%d\nM:%d\nB:%d\n",T,M,B);
  return 0;
}
