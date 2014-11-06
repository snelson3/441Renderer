#include <iostream>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <math.h>
#include "reader1F.cxx"
#include "shading.cxx"
#include "matrix.cxx"

using std::cerr;
using std::endl;

bool k = false;

double ceil441(double f)
{
    return ceil(f-0.00001);
}

double floor441(double f)
{
    return floor(f+0.00001);
}

vtkImageData *
NewImage(int width, int height)
{
    vtkImageData *img = vtkImageData::New();
    img->SetDimensions(width, height, 1);
    img->AllocateScalars(VTK_UNSIGNED_CHAR, 3);

    return img;
}

void
WriteImage(vtkImageData *img, const char *filename)
{
   std::string full_filename = filename;
   full_filename += ".png";
   vtkPNGWriter *writer = vtkPNGWriter::New();
   writer->SetInputData(img);
   writer->SetFileName(full_filename.c_str());
   writer->Write();
   writer->Delete();
}



double dotProduct(double *A, double *B)
{
  return A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
}

void debugOt(Point v1, Point v2, Point v3, std::vector<Triangle>::iterator i)
{
  cerr<<"***************"<<endl;
  cerr<<"v1("<<v1.X<<","<<v1.Y<<","<<v1.Z<<")"<<endl;
  cerr<<"v2("<<v2.X<<","<<v2.Y<<","<<v2.Z<<")"<<endl;
  cerr<<"v3("<<v3.X<<","<<v3.Y<<","<<v3.Z<<")"<<endl;
  cerr<<"**************"<<endl;
}


double findx (Point v1, Point v2, double E)
{
  return ((E-v1.Y)*(v2.X-v1.X))/(v2.Y-v1.Y)+v1.X;
}

double interpolate (double v1c, double v2c, double v1x, double v2x, double F)
{
  return v1c+((F-v1x)/(v2x-v1x))*(v2c-v1c);
}

class Screen
{
  public:
      unsigned char   *buffer;
      double *zbuffer;
      int width, height;
      LightingParameters lp;

      void setRed(double r, int x, int y)
      {
        int loc = y * (width)+x;
        buffer[loc*3] = (unsigned char)ceil441(fmin(1,r)*255.0);
      }

      void setGreen(double g, int x, int y)
      {
        int loc = y * (width) + x;
        buffer[(loc*3)+1] = (unsigned char)ceil441(fmin(1,g)*255.0);
      }
      void setBlue(double u, int x, int y){
        int loc = y * (width) + x;
        buffer[(loc*3)+2] = (unsigned char)ceil441(fmin(1,u)*255.0);
      }

      double getZ(int x, int y){
        int loc = y * (width) + x;
        return zbuffer[loc];
      }

      void setZ(int x, int y, double z){
        int loc = y *(width) + x;
        zbuffer[loc] = z;
      }

      double calculateDiffuse(double *normal)
      {
        double diff = fabs(dotProduct(lp.lightDir,normal));
        return diff;
      }

      double calculateSpecular(double *normal)
      {
        //View direction is going to be [0,0,-1];
        double V[3];
        V[0] = 0;
        V[1] = 0;
        V[2] = -1;
        double R[3];
        R[0] = 2 * dotProduct(lp.lightDir,normal)*normal[0]-lp.lightDir[0];
        R[1] = 2 * dotProduct(lp.lightDir,normal)*normal[1]-lp.lightDir[1];
        R[2] = 2 * dotProduct(lp.lightDir,normal)*normal[2]-lp.lightDir[2];
        return dotProduct(R,V);
      }

      double calculateShading(double x, double y, double *normal)
      {
        double diffuse = calculateDiffuse(normal);
        double specular = fmax(0,pow(calculateSpecular(normal),lp.alpha));
        double shading = lp.Ka + lp.Kd*diffuse + lp.Ks*specular;
        return shading;
      }

      void setColor(double x, double y, Point v1, Point v2, Point v3)
      {
        if (x == this->width)
            return;
        double v5x = findx(v1,v2,y);
        double v6x = findx(v2,v3,y);

        double v5z = interpolate(v1.Z,v2.Z,v1.Y,v2.Y,y);
        double v6z = interpolate(v2.Z,v3.Z,v2.Y,v3.Y,y);
        double v4z = interpolate(v5z,v6z,v5x,v6x,x);

         if ((v4z <= getZ(x,y) ))
            return;

        setZ(x,y,v4z);

        double v5r = interpolate(v1.color[0],v2.color[0],v1.Y,v2.Y,y);
        double v6r = interpolate(v2.color[0],v3.color[0],v2.Y,v3.Y,y);
        double v4r = interpolate(v5r,v6r,v5x,v6x,x);

        double v5u = interpolate(v1.color[2],v2.color[2],v1.Y,v2.Y,y);
        double v6u = interpolate(v2.color[2],v3.color[2],v2.Y,v3.Y,y);
        double v4u = interpolate(v5u,v6u,v5x,v6x,x);

        double v5g = interpolate(v1.color[1],v2.color[1],v1.Y,v2.Y,y);
        double v6g = interpolate(v2.color[1],v3.color[1],v2.Y,v3.Y,y);
        double v4g = interpolate(v5g,v6g,v5x,v6x,x);

        double v5n[3];
        double v6n[3];
        double v4n[3];

        v5n[0] = interpolate(v1.normal[0],v2.normal[0],v1.Y,v2.Y,y);
        v5n[1] = interpolate(v1.normal[1],v2.normal[1],v1.Y,v2.Y,y);
        v5n[2] = interpolate(v1.normal[2],v2.normal[2],v1.Y,v2.Y,y);

        v6n[0] = interpolate(v2.normal[0],v3.normal[0],v2.Y,v3.Y,y);
        v6n[1] = interpolate(v2.normal[1],v3.normal[1],v2.Y,v3.Y,y);
        v6n[2] = interpolate(v2.normal[2],v3.normal[2],v2.Y,v3.Y,y);

        v4n[0] = interpolate(v5n[0],v6n[0],v5x,v6x,x);
        v4n[1] = interpolate(v5n[1],v6n[1],v5x,v6x,x);
        v4n[2] = interpolate(v5n[2],v6n[2],v5x,v6x,x);

        this->setRed(v4r*this->calculateShading(x,y,v4n), x, y);
        this->setBlue(v4u*this->calculateShading(x,y,v4n), x, y);
        this->setGreen(v4g*this->calculateShading(x,y,v4n), x, y);
        }
};

void writeFlatBottom(Point v1, Point v2, Point v3, Screen screen)
{
    double leftslope, rightslope, li, ri, leftend, rightend;

    /********FIND THE LEFT AND RIGHT SLOPE AND INTERCEPT***********/
    if (v1.X==v2.X)
      leftslope = 0;//right triangle
    else
      leftslope = (v2.Y-v1.Y)/(v2.X-v1.X);
    li = v1.Y-leftslope*v1.X;

    if (v2.X == v3.X)
      rightslope = 0; //right triangle
    else
      rightslope = (v2.Y - v3.Y)/(v2.X - v3.X);
    ri = v3.Y - v3.X*rightslope;

    /****Loop through every hor line of the triangle***/

    for (int y = ceil441(v1.Y); y <= floor441(v2.Y); y++)
    {
      if (y >= 0 && y < screen.height)
      {
        //now I need to find the left endpoint and the right endpoint
        if (leftslope == 0) leftend = v1.X;
        else leftend = (y - li) / leftslope;
        if (rightslope == 0) rightend = v3.X;
        else rightend = (y - ri) / rightslope;

        if (ceil441(leftend) == floor(rightend) && (leftend >= 0 && leftend < screen.width))
            screen.setColor(ceil(leftend),y,v1,v2,v3);

        for (int x = ceil441(leftend); x <= floor441(rightend); x++)
          if (x >= 0 && x< screen.width)
            screen.setColor(x,y,v1,v2,v3);
      }
    }
}

void writeFlatTop(Point v1, Point v2, Point v3, Screen screen)
{
    double leftslope, rightslope, li, ri,leftend,rightend;

    if (v1.X==v2.X)
      leftslope = 0;//right triangle
    else
      leftslope = (v2.Y-v1.Y)/(v2.X-v1.X);
    li = v1.Y-leftslope*v1.X;

    if (v2.X == v3.X)
      rightslope = 0; //right triangle
    else
      rightslope = (v2.Y - v3.Y)/(v2.X - v3.X);
    ri = v3.Y - v3.X * rightslope;

    /****Loop through every hor line of the Triangle***/

    for (int y = floor441(v1.Y); y >= ceil441(v2.Y); y--)
    {
      if (y >= 0 && y < screen.height)
      {
      //now I need to find the left endpoint and the right endpoint
        if (leftslope == 0) leftend = v1.X;
        else leftend = (y - li) / leftslope;
        if (rightslope == 0) rightend = v3.X;
        else rightend = (y - ri) / rightslope;

        if ( ceil441(leftend) == floor441(rightend) && ( (leftend >= 0) && (leftend< screen.width) ) )
          screen.setColor(ceil441(leftend),y,v1,v2,v3);

        for (int x = ceil441(leftend); x <= floor441(rightend); x++)
          if (x >= 0 && x< screen.width)
            screen.setColor(x,y,v1,v2,v3);
      }
    }
}


void depositTriangle(std::vector<Triangle>::iterator i, Screen screen)
{
  Point v1,v2,v3,v4;
  double slope, intercept;

  /*Step 1*/

  //3 Cases, when 1 and 2 are flat bottom and when 1 and 3 are flat bottom and when 2 and 3 are flat bottom
  if (i->getY1() == i->getY2())
  {
    if (i->getX1() < i->getX2())
      i->setV(&v1,&v3,&v2);
    else
      i->setV(&v3,&v1,&v2);

    if (v2.Y < v1.Y)
      writeFlatTop(v1,v2,v3,screen);
    else
      writeFlatBottom(v1,v2,v3,screen);
  }
  else if (i->getY2() == i->getY3())
  {
    if (i->getX2() < i->getX3())
      i->setV(&v2,&v1,&v3);
    else
      i->setV(&v2,&v3,&v1);
    
    if (v2.Y < v1.Y)
      writeFlatTop(v1,v2,v3,screen);
    else
      writeFlatBottom(v1,v2,v3,screen);
  }
  else if (i->getY1() == i->getY3())
  {
    if (i->getX1() < i->getX3())
      i->setV(&v1,&v2,&v3);
    else
      i->setV(&v3,&v2,&v1);

    if (v2.Y < v1.Y)
      writeFlatTop(v1,v2,v3,screen);
    else
      writeFlatBottom(v1,v2,v3,screen);
  }
  else
  {
    //must split into a flat bottom and a flat top triangle

    /*Step 2*/


    i->sortVertices(&v1,&v2,&v3);

    /*Step 3*/

    if (v1.X == v3.X) slope = 0; //right triangle
    else slope = (v1.Y - v3.Y)/(v1.X - v3.X);
    intercept = v3.Y - v3.X*slope;

    /*Step 4*/

    v4.Y = v2.Y;

    if (slope == 0) v4.X = v3.X;
    else v4.X = (v4.Y - intercept) / slope;

    /*Step 4.5*/
    
    v4.Z = interpolate(v1.Z,v3.Z,v1.Y,v3.Y,v4.Y);

    v4.color[0] = interpolate(v1.color[0],v3.color[0],v1.Y,v3.Y,v4.Y);
    v4.color[1] = interpolate(v1.color[1],v3.color[1],v1.Y,v3.Y,v4.Y);
    v4.color[2] = interpolate(v1.color[2],v3.color[2],v1.Y,v3.Y,v4.Y);

    v4.normal[0] = interpolate(v1.normal[0],v3.normal[0],v1.Y,v3.Y,v4.Y);
    v4.normal[1] = interpolate(v1.normal[1],v3.normal[1],v1.Y,v3.Y,v4.Y);
    v4.normal[2] = interpolate(v1.normal[2],v3.normal[2],v1.Y,v3.Y,v4.Y);          
    
    /*Step 5*/

    if (v2.X < v4.X)
      writeFlatTop(v2,v3,v4,screen);
    else
      writeFlatTop(v4,v3,v2,screen);

    /*Step 6*/

    if (v2.X<v4.X)
      writeFlatBottom(v2,v1,v4,screen);
    else 
      writeFlatBottom(v4,v1,v2,screen);
  }
}



class Camera
{
  public:
    double          near, far;
    double          angle;
    double          position[3];
    double          focus[3];
    double          up[3];

    Matrix          ViewTransform(void) {Matrix vt; return vt;};
    Matrix          CameraTransform(void) {Matrix ct; return ct;};
    Matrix          DeviceTransform(void) {Matrix dt; return dt;};
};


double SineParameterize(int curFrame, int nFrames, int ramp)
{
    int nNonRamp = nFrames-2*ramp;
    double height = 1./(nNonRamp + 4*ramp/M_PI);
    if (curFrame < ramp)
    {
        double factor = 2*height*ramp/M_PI;
        double eval = cos(M_PI/2*((double)curFrame)/ramp);
        return (1.-eval)*factor;
    }
    else if (curFrame > nFrames-ramp)
    {
        int amount_left = nFrames-curFrame;
        double factor = 2*height*ramp/M_PI;
        double eval =cos(M_PI/2*((double)amount_left/ramp));
        return 1. - (1-eval)*factor;
    }
    double amount_in_quad = ((double)curFrame-ramp);
    double quad_part = amount_in_quad*height;
    double curve_part = height*(2*ramp)/M_PI;
    return quad_part+curve_part;
}

Camera
GetCamera(int frame, int nframes)
{
    double t = SineParameterize(frame, nframes, nframes/10);
    Camera c;
    c.near = 5;
    c.far = 200;
    c.angle = M_PI/6;
    c.position[0] = 40*sin(2*M_PI*t);
    c.position[1] = 40*cos(2*M_PI*t);
    c.position[2] = 40;
    c.focus[0] = 0;
    c.focus[1] = 0;
    c.focus[2] = 0;
    c.up[0] = 0;
    c.up[1] = 1;
    c.up[2] = 0;
    return c;
}






Triangle toDeviceSpace(std::vector<Triangle>::iterator i)
  {
    Triangle t;
    return t;
  }

void makeImage(int curr, int total, std::vector<Triangle> t)
{
  //set up the screen
vtkImageData *image = NewImage(1000, 1000);
   unsigned char *buffer =
     (unsigned char *) image->GetScalarPointer(0,0,0);

   int npixels = 1000*1000;
    double *zbuffer = new double[npixels];
    for (int i = 0; i < npixels*3 ; i++)
        buffer[i] = 0;
    for (int i = 0; i < npixels; i++)
        zbuffer[i] = 1; //I was told depth is opposite in this part

   Screen screen;
   screen.buffer = buffer;
   screen.zbuffer = zbuffer;
   screen.width = 1000;
   screen.height = 1000;

  //Now to get the image, put each triangle in device space, and render it on the screen (before writing the png)

  Camera c = GetCamera(curr,total);
  //I probably need to do something to get the right Matrix here to send into toDeviceSpace

  for(std::vector<Triangle>::iterator i = t.begin(); i != t.end(); ++i)
  {
    Triangle t = toDeviceSpace(i);
    depositTriangle(i,screen);
  }

  //This breakfs with big amounts of frames but thats ok for now/this project
  char framenum[6];
  if (curr < 10) sprintf(framenum, "fr000%d",curr);
  else if (curr < 100) sprintf(framenum, "fr00%d",curr);
  else if (curr < 1000) sprintf(framenum, "fr0%d",curr);
  else  sprintf(framenum, "fr%d",curr);
  WriteImage(image, framenum);
}

int main()
{
  std::vector<Triangle> triangles = GetTriangles();
  makeImage(0,1000,triangles);
  makeImage(250,1000,triangles);
  makeImage(500,1000,triangles);
  makeImage(750,1000,triangles);
   // vtkImageData *image = NewImage(1000, 1000);
   // unsigned char *buffer =
   //   (unsigned char *) image->GetScalarPointer(0,0,0);

   // int npixels = 1000*1000;
   //  double *zbuffer = new double[npixels];
   //  for (int i = 0; i < npixels*3 ; i++)
   //      buffer[i] = 0;
   //  for (int i = 0; i < npixels; i++)
   //      zbuffer[i] = -1; //valid depth values go from -1 (back) to 0 (front)

   // Screen screen;
   // screen.buffer = buffer;
   // screen.zbuffer = zbuffer;
   // screen.width = 1000;
   // screen.height = 1000;

   // for(std::vector<Triangle>::iterator i = triangles.begin(); i != triangles.end(); ++i)
   //   depositTriangle(i,screen);
   
   // WriteImage(image, "allTriangles");
}
