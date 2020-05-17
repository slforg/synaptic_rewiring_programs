#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

#include<X11/Xlib.h>
#include<X11/Xutil.h>
#include<unistd.h>

#define XSIZE 48
#define YSIZE 48
#define NM     8
#define DIRECTION 12
#define GRADE    10
#define COOR_MAX 120
#define WINDOW_X 100
#define WINDOW_Y 100
#define DOT_SIZE 5
#define WINDOW_WIDTH XSIZE*DOT_SIZE
#define WINDOW_HEIGHT YSIZE*DOT_SIZE
#define WINDOW_F 5 
#define clip 0.85

//Color Return Function ex GetColor(dis,"yellow");
unsigned long GetColor(Display* dis, char* color_name){
  Colormap cmap;
  XColor near_color,true_color;

  cmap=DefaultColormap(dis,0);
  XAllocNamedColor(dis,cmap,color_name,&near_color,&true_color);
  return (near_color.pixel);
}

//Dot Color
void DotRGB(Display *dis,Window win,GC gc,int x,int y,unsigned char r,unsigned char g,unsigned char b){
  XSetForeground(dis,gc,(r<<16)|(g<<8)|(b)<<0);
  XFillRectangle(dis,win,gc,x,y,DOT_SIZE,DOT_SIZE);
}

main(){
  int x,y,x0,y0,dx,dy,xd,yd,i0,j0,k,l,n,i,j;
  float r, mgn;
  int sx,sy; //select x&y point
  int dir_i;
  char name[10];
  int color_Set[256][256];

  float dir[XSIZE][YSIZE][NM];
  Display *dis;
  Window rw,w;
  GC gc;
  XSetWindowAttributes att;

  typedef struct{
    unsigned int r,g,b;
  }COLOR;

  COLOR rgb[DIRECTION][GRADE];

  FILE *fdir;

  printf("input file name : ");
  scanf("%s",name);
  getchar();

  if((fdir = fopen(name,"rb"))==NULL){
    printf("Can not open output file.\n");
    exit(-1);
  }
  for(j0=0;j0<YSIZE;j0++){
    for(i0=0;i0<YSIZE;i0++){
     fscanf(fdir,"%f %f %f %f %f %f %f %f\n"
			       ,&dir[i0][j0][0],&dir[i0][j0][1]
			       ,&dir[i0][j0][2],&dir[i0][j0][3]
			       ,&dir[i0][j0][4],&dir[i0][j0][5]
			       ,&dir[i0][j0][6],&dir[i0][j0][7]);
    }
  }

  dis = XOpenDisplay(NULL);
  rw =XDefaultRootWindow(dis);
 
  w=XCreateSimpleWindow(dis,
			rw,
			WINDOW_X,
			WINDOW_Y,
			WINDOW_WIDTH,
			WINDOW_HEIGHT,
			WINDOW_F,
			WhitePixel(dis,0),
			BlackPixel(dis,0));
  XStoreName(dis,w,"Direction_Map"); //window name
  att.backing_store=Always;
  XChangeWindowAttributes(dis,w,CWBackingStore,&att);
  XMapWindow(dis,w); //Window Map
  gc = XCreateGC(dis,rw,0,0);
  
  for(k=0;k<GRADE;k++){
    r=(float)k/10.0;
    l=(10-k)*6;
    //0
    rgb[0][k].r=(int)(0*r);
    rgb[0][k].g=(int)(0*r);
    rgb[0][k].b=(int)(255-255*r);
    //30
    rgb[1][k].r=(int)(0*r);
    rgb[1][k].g=(int)(255-255*r);
    rgb[1][k].b=(int)(255-255*r);
    //60
    rgb[2][k].r=(int)(0*r);
    rgb[2][k].g=(int)(255-255*r);
    rgb[2][k].b=(int)(0*r);
    //90
    rgb[3][k].r=(int)(255-255*r);
    rgb[3][k].g=(int)(255-255*r);
    rgb[3][k].b=(int)(0*r);
    //135
    rgb[4][k].r=(int)(255-255*r);
    rgb[4][k].g=(int)(0*r);
    rgb[4][k].b=(int)(0*r);
    //150
    rgb[5][k].r=(int)(255-255*r);
    rgb[5][k].g=(int)(0*r);
    rgb[5][k].b=(int)(255-255*r);
  }
  
  for(y0=0;y0<YSIZE;y0++){
    for(x0=0;x0<XSIZE;x0++){
      dir_i=(int)(dir[x0][y0][5]);
      if(165<=dir_i || dir_i< 15){i=0;}
      if( 15<=dir_i && dir_i< 45){i=1;}
      if( 45<=dir_i && dir_i< 75){i=2;}
      if( 75<=dir_i && dir_i<105){i=3;}
      if(105<=dir_i && dir_i<135){i=4;}
      if(135<=dir_i && dir_i<165){i=5;}
      mgn=dir[x0][y0][6]/clip;
      if(1.0<=mgn){mgn=1.0;}
      j=9-(int)(mgn*9);
      //      printf("%f \n", mgn);
      DotRGB(dis,w,gc,x0*DOT_SIZE,WINDOW_HEIGHT-(y0+1)*DOT_SIZE,rgb[i][j].r,rgb[i][j].g,rgb[i][j].b);
    }
  }

  

  XFlush(dis);
  system("import -window Direction_Map map.bmp");
  
  printf("Hit Return Key\n");
  getchar();
  XFreeGC(dis,gc);
  XDestroyWindow(dis,w);
  XFlush(dis);
  XCloseDisplay(dis);
}

