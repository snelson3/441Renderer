#include <iostream>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
//#include "get_triangles.cxx"
#include <vtkPointData.h>
#include <math.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>


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

struct Point
{
  double X;
  double Y;
  double Z;
};


class Triangle
{
  public:
      double         X[3];
      double         Y[3];
      double         Z[3];
      double        colors[3][3];

      Point *PT1;
      Point *PT2;
      Point *PT3;

      void setPT1(double x, double y, double z){X[0]=x; Y[0]=y; Z[0] = z;};
      void setPT2(double x, double y, double z){X[1]=x; Y[1]=y; Z[1] = z;};
      void setPT3(double x, double y, double z){X[2]=x; Y[2]=y; Z[2] = z;};

      Point getPT(int t)
      {
        Point tr;
        tr.X = X[t];
        tr.Y = Y[t];
        tr.Z = Z[t];
        return tr;
      }

      double getX1(){return X[0];};
      double getX2(){return X[1];};
      double getX3(){return X[2];};

      double getY1(){return Y[0];};
      double getY2(){return Y[1];};
      double getY3(){return Y[2];};

      double getZ1(){return Z[0];};
      double getZ2(){return Z[1];};
      double getZ3(){return Z[2];};

      //This is not the correct places to put the colors I think
        //I'm not even sure why the pixels are a 2d array now.
      void setRed(double c,int i){colors[i][0] = c;};
      void setBlue(double c, int i){colors[i][1] = c;};
      void setGreen(double c, int i){colors[i][2] = c;};

      //these need to be changed too.
      double getRed(int i){return colors[i][0];};
      double getBlue(int i){return colors[i][1];};
      double getGreen(int i){return colors[i][2];};
 
      void debugOut()
      {
            std::cerr<<"********************************"<<endl;
            std::cerr<<"V1("<<this->getX1()<<","<<this->getY1()<<","<<this->getZ1()<<")"<<endl;
            std::cerr<<"V2("<<this->getX2()<<","<<this->getY2()<<","<<this->getZ2()<<")"<<endl;
            std::cerr<<"V3("<<this->getX3()<<","<<this->getY3()<<","<<this->getZ3()<<")"<<endl;
            //std::cerr<<"C("<<this->getRed()<<","<<this->getBlue()<<","<<this->getGreen()<<")"<<endl;
      }

      void setV(Point *v1, Point *v2, Point *v3)
      {
        *v1 = getPT(0);
        *v2 = getPT(1);
        *v3 = getPT(2);
        PT1 = v1;
        PT2 = v2;
        PT3 = v3;//I think its because
      }

      void sortVertices(Point v1, Point v2, Point v3)
      {
        double vt1;
        double vt3;
        //could I use a bunch of for loops to reduce the logic here?
        vt1 = fmax(Y[0],fmax(Y[1],Y[2]));
        vt3 = fmin(Y[0],fmin(Y[1],Y[2]));
        if (vt1 == Y[0] && vt3 == Y[1])
        {
          setV(&v1,&v3,&v2);
        }
        else if (vt1 == Y[0] && vt3 == Y[2])
        {
          setV(&v1,&v2,&v3);

        }
        else if (vt1 == Y[1] && vt3 == Y[0])
        {
          setV(&v3,&v1,&v2);
        }
        else if (vt1 == Y[1] && vt3 == Y[2])
        {
          setV(&v2,&v1,&v3);
        }
        else if (vt1 == Y[2] && vt3 == Y[0])
        {
          setV(&v3,&v2,&v1);
        }
        else if (vt1 == Y[2] && vt3 == Y[1])
        {
          setV(&v2,&v3,&v1);
        }
        else
          std::cerr<<"SOMETHINGS BROkeN"<<endl;
      }

  // would some methods for transforming the triangle in place be helpful?
};

void debugOt(Point v1, Point v2, Point v3, std::vector<Triangle>::iterator i)
{
  cerr<<"***************"<<endl;
  cerr<<"v1("<<v1.X<<","<<v1.Y<<","<<v1.Z<<")"<<endl;
  cerr<<"v2("<<v2.X<<","<<v2.Y<<","<<v2.Z<<")"<<endl;
  cerr<<"v3("<<v3.X<<","<<v3.Y<<","<<v3.Z<<")"<<endl;
  cerr<<"**************"<<endl;
  cerr<<"TRIANGLE OF HELL";
  cerr<<endl<<"v1("<<i->PT1.X<<","<<i->PT1.Y<<","<<i->PT1.Z<<")"<<endl;
  cerr<<"v2("<<i->PT2.X<<","<<i->PT2.Y<<","<<i->PT2.Z<<")"<<endl;
  cerr<<"v3("<<i->PT3.X<<","<<i->PT3.Y<<","<<i->PT3.Z<<")"<<endl;
}


double findx (Point v1, Point v2, double E)
{
  return ((E-v1.Y)*(v2.X-v1.X))/(v2.Y-v1.Y)+v1.X;
}

double findc (double v1c, double v2c, double v1x, double v2x, double F)
{
  return v1c+((F-v1x)/(v2x-v1x))*(v2c-v1c);
}

class Screen
{
  public:
      unsigned char   *buffer;
      int width, height;

      //assume 9 pixels, so setr[3][3] would be spot 24

      void setRed(double r, int x, int y){
        //std::cerr<<"  SETRED "<<x<<endl;
        int loc = y * (width)+x;
        buffer[loc*3] = r*255;
        //buffer[(y*width)+x] = r;
      }
      void setBlue(double u, int x, int y){
        int loc = y * (width) + x;
        buffer[(loc*3)+1] = u*255;
        //buffer[(y*width)+x+1] = u;
      }
      void setGreen(double g, int x, int y){
        int loc = y * (width) + x;
        //buffer[(y*width)+x+2]=g;
        buffer[(loc*3)+2] = g*255;
        //cerr<<x<<y<<"\n";
      }

      //bool l = false;

      void setColor(std::vector<Triangle>::iterator i, double x, double y)//int x, int y)
      {
        //This is probably where I want to do the interpolation, but the scanline maybe
        // cerr<<" ";
        //cerr<<i->getRed(0,0)<<i->getBlue(0,0)<<i->getGreen(0,0);
        if (x == this->width)
          {
            return;
          }
        else
        {
          // cerr<<"SETTING COLOR"<<endl;
        //I think interpolating is basically slope

          //Interpolate red
        double v5x = findx(i->PT1,i->PT2,y);
        // cerr<<"V5x "<<v5x<<endl;
        double v6x = findx(i->PT2,i->PT3,y);
        // cerr<<"V6x "<<v6x<<endl;

        double v5r = findc(i->getRed(0),i->getRed(1),i->PT1.X,i->PT2.X,v5x);

        double v6r = findc(i->getRed(1),i->getRed(2),i->PT2.X,i->PT3.X,v6x);
        double v4r = findc(v5r,v6r,v5x,v6x,x);

        // cerr<<"V5r"<<v5r<<endl;
        // cerr<<" A "<<i->getRed(0)<<endl;
        // cerr<<" B "<<i->getRed(1)<<endl;
        // cerr<<" C "<<i->PT1.X<<endl;
        // cerr<<" D "<<i->PT2.X<<endl;
        // cerr<<" F "<<v5x<<endl;
        // cerr<<"V6r"<<v6r<<endl;


        double v5u = findc(i->getBlue(0),i->getBlue(1),i->PT1.X,i->PT2.X,v5x);
        double v6u = findc(i->getBlue(1),i->getBlue(2),i->PT2.X,i->PT3.X,v6x);
        double v4u = findc(v5u,v6u,v5x,v6x,x);

        double v5g = findc(i->getGreen(0),i->getGreen(1),i->PT1.X,i->PT2.X,v5x);
        double v6g = findc(i->getGreen(1),i->getGreen(2),i->PT2.X,i->PT3.X,v6x);
        double v4g = findc(v5g,v6g,v5x,v6x,x);

         this->setRed(v4r, x, y);
         this->setBlue(v4u, x, y);
         this->setGreen(v4g, x, y);

        //this->setRed(111,x,y);
        //this->setBlue(111,x,y);
        //this->setGreen(111,x,y);

        // cerr<<"V4R"<<v4r;
        }
      }


};


std::vector<Triangle>
GetTriangles(void)
{
    vtkPolyDataReader *rdr = vtkPolyDataReader::New();
    rdr->SetFileName("proj1d_geometry.vtk");
    cerr << "Reading" << endl;
    rdr->Update();
    cerr << "Done reading" << endl;
    if (rdr->GetOutput()->GetNumberOfCells() == 0)
    {
        cerr << "Unable to open file!!" << endl;
        exit(EXIT_FAILURE);
    }
    vtkPolyData *pd = rdr->GetOutput();
    int numTris = pd->GetNumberOfCells();
    vtkPoints *pts = pd->GetPoints();
    vtkCellArray *cells = pd->GetPolys();
    vtkFloatArray *var = (vtkFloatArray *) pd->GetPointData()->GetArray("hardyglobal");
    float *color_ptr = var->GetPointer(0);
    std::vector<Triangle> tris(numTris);
    vtkIdType npts;
    vtkIdType *ptIds;
    int idx;
    for (idx = 0, cells->InitTraversal() ; cells->GetNextCell(npts, ptIds) ; idx++)
    {
        if (npts != 3)
        {
            cerr << "Non-triangles!! ???" << endl;
            exit(EXIT_FAILURE);
        }
        tris[idx].X[0] = pts->GetPoint(ptIds[0])[0];
        tris[idx].X[1] = pts->GetPoint(ptIds[1])[0];
        tris[idx].X[2] = pts->GetPoint(ptIds[2])[0];
        tris[idx].Y[0] = pts->GetPoint(ptIds[0])[1];
        tris[idx].Y[1] = pts->GetPoint(ptIds[1])[1];
        tris[idx].Y[2] = pts->GetPoint(ptIds[2])[1];
        tris[idx].Z[0] = pts->GetPoint(ptIds[0])[2];
        tris[idx].Z[1] = pts->GetPoint(ptIds[1])[2];
        tris[idx].Z[2] = pts->GetPoint(ptIds[2])[2];
        // 1->2 interpolate between light blue, dark blue
        // 2->2.5 interpolate between dark blue, cyan
        // 2.5->3 interpolate between cyan, green
        // 3->3.5 interpolate between green, yellow
        // 3.5->4 interpolate between yellow, orange
        // 4->5 interpolate between orange, brick
        // 5->6 interpolate between brick, salmon
        double mins[7] = { 1, 2, 2.5, 3, 3.5, 4, 5 };
        double maxs[7] = { 2, 2.5, 3, 3.5, 4, 5, 6 };
        unsigned char RGB[8][3] = { { 71, 71, 219 }, 
                                    { 0, 0, 91 },
                                    { 0, 255, 255 },
                                    { 0, 128, 0 },
                                    { 255, 255, 0 },
                                    { 255, 96, 0 },
                                    { 107, 0, 0 },
                                    { 224, 76, 76 } 
                                  };
        for (int j = 0 ; j < 3 ; j++)
        {
            float val = color_ptr[ptIds[j]];
            int r;
            for (r = 0 ; r < 7 ; r++)
            {
                if (mins[r] <= val && val < maxs[r])
                    break;
            }
            if (r == 7)
            {
                cerr << "Could not interpolate color for " << val << endl;
                exit(EXIT_FAILURE);
            }
            double proportion = (val-mins[r]) / (maxs[r]-mins[r]);
            tris[idx].colors[j][0] = (RGB[r][0]+proportion*(RGB[r+1][0]-RGB[r][0]))/255.0;
            tris[idx].colors[j][1] = (RGB[r][1]+proportion*(RGB[r+1][1]-RGB[r][1]))/255.0;
            tris[idx].colors[j][2] = (RGB[r][2]+proportion*(RGB[r+1][2]-RGB[r][2]))/255.0;
        }
    }

    return tris;
}



void writeFlatBottom(Point v1, Point v2, Point v3, Screen screen,
              std::vector<Triangle>::iterator i)
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
            screen.setColor(i,ceil(leftend),y);

        for (int x = ceil441(leftend); x <= floor441(rightend); x++)
          if (x >= 0 && x< screen.width)
            screen.setColor(i,x,y);
      }
    }
}

void writeFlatTop(Point v1, Point v2, Point v3, Screen screen,
                    std::vector<Triangle>::iterator i)
{
    double leftslope, rightslope, li, ri,leftend,rightend;
    //debugOt(v1,v2,v3,i);
    /***Find the Left and Right Slope and Intercept***/

    //get rid of this
    //i->setBlue(255);
    //i->setGreen(0);
    //i->setRed(0);

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
          screen.setColor(i,ceil441(leftend),y);

        for (int x = ceil441(leftend); x <= floor441(rightend); x++)
          if (x >= 0 && x< screen.width)
            screen.setColor(i,x,y);
      }
    }
}


      void depositTriangle(std::vector<Triangle>::iterator i, Screen screen)
      {
        Point v1,v2,v3,v4;
        double slope, intercept;

        /*Step 1*/

    /*****I might need to do a continue if three vertices are in a line but maybe not*******/


        //3 Cases, when 1 and 2 are flat bottom and when 1 and 3 are flat bottom and when 2 and 3 are flat bottom
        if (i->getY1() == i->getY2())
        {
          //I think these can fuck up given a triangle like (0,0,1)(0,0,2)(5,2,0)
          //1 and 2 are the same Y
          if (i->getX1() < i->getX2())
          {
            i->setV(&v1,&v3,&v2);
          }
          else
          {
            i->setV(&v2,&v3,&v1);
          }
          if (v2.Y < v1.Y)
            writeFlatTop(v1,v2,v3,screen,i);
          else
            writeFlatBottom(v1,v2,v3,screen,i);
        }
        else if (i->getY2() == i->getY3())
        {
          //2 and 3 are the same Y
          if (i->getX2() < i->getX3())
          {
            i->setV(&v2,&v1,&v3);
          }
          else
          {
            i->setV(&v3,&v1,&v2);
          }
          if (v2.Y < v1.Y)
            writeFlatTop(v1,v2,v3,screen,i);
          else
            writeFlatBottom(v1,v2,v3,screen,i);
        }
        else if (i->getY1() == i->getY3())
        {    
          //1 and 3 are the same Y,
          if (i->getX1() < i->getX3())
          {
            i->setV(&v1,&v2,&v3);
          }
          else
          {
            i->setV(&v3,&v2,&v1);
          }
          if (v2.Y < v1.Y)
            writeFlatTop(v1,v2,v3,screen,i);
          else
          {
            writeFlatBottom(v1,v2,v3,screen,i);
          }
        }
        else
        {  
          //must split into a flat bottom and a flat top triangle

          /*Step 2*/

          i->sortVertices(v1,v2,v3);

          /*Step 3*/

          if (v1.X == v3.X) slope = 0; //right triangle
          else slope = (v1.Y - v3.Y)/(v1.X - v3.X);
          intercept = v3.Y - v3.X*slope;

          /*Step 4*/

          v4.Y = v2.Y;

          if (slope == 0) v4.X = v3.X;
          else v4.X = (v4.Y - intercept) / slope;

          /*Step 4.5*/
          //Gotta do some math to determine the Z for this point
          v4.Z = 0;

          /*Step 5*/

          if (v2.X < v4.X)
            writeFlatTop(v2,v3,v4,screen,i);
          else
            writeFlatTop(v4,v3,v2,screen,i);

          /*Step 6*/
        
          if (v2.X<v4.X)
            writeFlatBottom(v2,v1,v4,screen,i);
          else 
            writeFlatBottom(v4,v1,v2,screen,i);

    }
      }

int main()
{
   vtkImageData *image = NewImage(1000, 1000);
   unsigned char *buffer =
     (unsigned char *) image->GetScalarPointer(0,0,0);
   int npixels = 1000*1000;
   for (int i = 0 ; i < npixels*3 ; i++)
       buffer[i] = 0;

   std::vector<Triangle> triangles = GetTriangles();

   Screen screen;
   screen.buffer = buffer;
   screen.width = 1000;
   screen.height = 1000;

   // YOUR CODE GOES HERE TO DEPOSIT TRIANGLES INTO PIXELS USING THE SCANLINE ALGORITHM
   int e = 0;
  for(std::vector<Triangle>::iterator i = triangles.begin(); i != triangles.end(); ++i) 
  {
    if (e == 1) continue;
    i->X[0] = 200; i->X[1] = 700; i->X[2] = 400;
    i->Y[0] = 500; i->Y[1] = 500; i->Y[2] = 300; //I think it's something to do with the y value
    i->Z[0] = 100; i->Z[1] = 100; i->Z[2] = 100;
    depositTriangle(i,screen);
    e++;
  }
  std::cerr<<e<<" ";


  std::cerr<<"Now IT's time To write";
   WriteImage(image, "allTriangles");
}
