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

class Triangle
{
  public:
      double         X[3];
      double         Y[3];
      double         Z[3];
      double        colors[3][3];

      void setPT1(double x, double y, double z){X[0]=x; Y[0]=y; Z[0] = z;};
      void setPT2(double x, double y, double z){X[1]=x; Y[1]=y; Z[1] = z;};
      void setPT3(double x, double y, double z){X[2]=x; Y[2]=y; Z[2] = z;};

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
      void setRed(double c,int i, int j){for (int i = 0; i<3;i++)colors[i][j] = c;};
      void setBlue(double c, int i, int j){for (int i=0; i<3;i++)colors[i][j] = c;};
      void setGreen(double c, int i, int j){for (int i=0;i<3;i++)colors[i][j] = c;};

      //these need to be changed too.
      unsigned char getRed(int i, int j){return colors[i][j];};
      unsigned char getBlue(int i, int j){return colors[i][j];};
      unsigned char getGreen(int i, int j){return colors[i][j];};

      void debugOut()
      {
            std::cerr<<"********************************"<<endl;
            std::cerr<<"V1("<<this->getX1()<<","<<this->getY1()<<","<<this->getZ1()<<")"<<endl;
            std::cerr<<"V2("<<this->getX2()<<","<<this->getY2()<<","<<this->getZ2()<<")"<<endl;
            std::cerr<<"V3("<<this->getX3()<<","<<this->getY3()<<","<<this->getZ3()<<")"<<endl;
            //std::cerr<<"C("<<this->getRed()<<","<<this->getBlue()<<","<<this->getGreen()<<")"<<endl;
      }

      void sortVertices(double v1[], double v2[], double v3[])
      {
        //could I use a bunch of for loops to reduce the logic here?
        v1[1] = fmax(Y[0],fmax(Y[1],Y[2]));
        v3[1] = fmin(Y[0],fmin(Y[1],Y[2]));
        if (v1[1] == Y[0] && v3[1] == Y[1])
        {
          v1[0] = X[0];
          v3[0] = X[1];
          v2[0] = X[2];
          v2[1] = Y[2];
          v1[2] = Z[0];
          v2[2] = Z[2];
          v3[2] = Z[1];
        }
        else if (v1[1] == Y[0] && v3[1] == Y[2])
        {
          v1[0] = X[0];
          v3[0] = X[2];
          v2[0] = X[1];
          v2[1] = Y[1];
          v1[2] = Z[0];
          v2[2] = Z[1];
          v3[2] = Z[2];

        }
        else if (v1[1] == Y[1] && v3[1] == Y[0])
        {
          v1[0] = X[1];
          v3[0] = X[0];
          v2[0] = X[2];
          v2[1] = Y[2];

          v1[2] = Z[1];
          v2[2] = Z[2];
          v3[2] = Z[0];
        }
        else if (v1[1] == Y[1] && v3[1] == Y[2])
        {
          v1[0] = X[1];
          v3[0] = X[2];
          v2[0] = X[0];
          v2[1] = Y[0];

          v1[2] = Z[1];
          v2[2] = Z[0];
          v3[2] = Z[2];
        }
        else if (v1[1] == Y[2] && v3[1] == Y[0])
        {
          v1[0] = X[2];
          v3[0] = X[0];
          v2[0] = X[1];
          v2[1] = Y[1];

          v1[2] = Z[2];
          v2[2] = Z[1];
          v3[2] = Z[0];
        }
        else if (v1[1] == Y[2] && v3[1] == Y[1])
        {
          v1[0] = X[2];
          v3[0] = X[1];
          v2[0] = X[0];
          v2[1] = Y[0];
          v1[2] = Z[2];
          v2[2] = Z[0];
          v3[2] = Z[1];
        }
        else
          std::cerr<<"SOMETHINGS BROkeN"<<endl;
      }






  // would some methods for transforming the triangle in place be helpful?
};

class Screen
{
  public:
      unsigned char   *buffer;
      int width, height;

      //assume 9 pixels, so setr[3][3] would be spot 24

      void setRed(double r, int x, int y){
        //std::cerr<<"  SETRED "<<x<<endl;
        int loc = y * (width)+x;
        buffer[loc*3] = r;
        //buffer[(y*width)+x] = r;
      }
      void setBlue(double u, int x, int y){
        int loc = y * (width) + x;
        buffer[(loc*3)+1] = u;
        //buffer[(y*width)+x+1] = u;
      }
      void setGreen(double g, int x, int y){
        int loc = y * (width) + x;
        //buffer[(y*width)+x+2]=g;
        buffer[(loc*3)+2] = g;
      }

bool l = false;

      void setColor(std::vector<Triangle>::iterator i, int x, int y)
      {
        //This is probably where I want to do the interpolation, but the scanline maybe
        if (x == this->width)
          {
            return;
          }
        else
        {
        this->setRed(i->getRed(0,0), x, y);
        this->setBlue(i->getBlue(0,0), x, y);
        this->setGreen(i->getGreen(0,0), x, y);
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



void writeFlatBottom(double v1[], double v2[], double v3[], Screen screen,
              std::vector<Triangle>::iterator i)
{
    double leftslope, rightslope, li, ri, leftend, rightend;

    /********FIND THE LEFT AND RIGHT SLOPE AND INTERCEPT***********/

    if (v1[0]==v2[0])
      leftslope = 0;//right triangle
    else
      leftslope = (v2[1]-v1[1])/(v2[0]-v1[0]);
    li = v1[1]-leftslope*v1[0];

    if (v2[0] == v3[0])
      rightslope = 0; //right triangle
    else
      rightslope = (v2[1] - v3[1])/(v2[0] - v3[0]);
    ri = v3[1] - v3[0]*rightslope;

    /****Loop through every hor line of the triangle***/

    for (int y = ceil441(v1[1]); y <= floor441(v2[1]); y++)
    {
      if (y >= 0 && y < screen.height)
      {
        //now I need to find the left endpoint and the right endpoint
        if (leftslope == 0) leftend = v1[0];
        else leftend = (y - li) / leftslope;
        if (rightslope == 0) rightend = v3[0];
        else rightend = (y - ri) / rightslope;

        if (ceil441(leftend) == floor(rightend) && (leftend >= 0 && leftend < screen.width))
            screen.setColor(i,ceil(leftend),y);

        for (int x = ceil441(leftend); x <= floor441(rightend); x++)
          if (x >= 0 && x< screen.width)
            screen.setColor(i,x,y);
      }
    }
}

void writeFlatTop(double v1[], double v2[], double v3[], Screen screen,
                    std::vector<Triangle>::iterator i)
{
    double leftslope, rightslope, li, ri,leftend,rightend;

    /***Find the Left and Right Slope and Intercept***/

    //get rid of this
    //i->setBlue(255);
    //i->setGreen(0);
    //i->setRed(0);

    if (v1[0]==v2[0])
      leftslope = 0;//right triangle
    else
      leftslope = (v2[1]-v1[1])/(v2[0]-v1[0]);
    li = v1[1]-leftslope*v1[0];

    if (v2[0] == v3[0])
      rightslope = 0; //right triangle
    else
      rightslope = (v2[1] - v3[1])/(v2[0] - v3[0]);
    ri = v3[1] - v3[0] * rightslope;

    /****Loop through every hor line of the Triangle***/

    for (int y = floor441(v1[1]); y >= ceil441(v2[1]); y--)
    {
      if (y >= 0 && y < screen.height)
      {
      //now I need to find the left endpoint and the right endpoint
        if (leftslope == 0) leftend = v1[0];
        else leftend = (y - li) / leftslope;
        if (rightslope == 0) rightend = v3[0];
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
        double v1 [3]; //v1 is bottom left vertex
        double v2 [3]; //v2 is top left vertex
        double v3 [3]; //v3 is bottom right vertex
        double v4 [3];
        double slope, intercept;

        /*Step 1*/

    /*****I might need to do a continue if three vertices are in a line but maybe not*******/



        //3 Cases, when 1 and 2 are flat bottom and when 1 and 3 are flat bottom and when 2 and 3 are flat bottom
        if (i->getY1() == i->getY2())
        {
          //I think these can fuck up given a triangle like (0,0,1)(0,0,2)(5,2,,0)
          //1 and 2 are the same Y
          if (i->getX1() < i->getX2())
          {
            v1[0] = i->getX1(); v1[1] = i->getY1(); v1[2] = i->getZ1();
            v3[0] = i->getX2(); v3[1] = i->getY2(); v1[2] 
            v2[0] = i->getX3();
          }
          else
          {
            v1[0] = i->getX2();
            v3[0] = i->getX1();
            v2[0] = i->getX3();
          }
          if (v2[1] < v1[1])
            writeFlatTop(v1,v2,v3,screen,i);
          else
            writeFlatBottom(v1,v2,v3,screen,i);
        }
        else if (i->getY2() == i->getY3())
        {
          //2 and 3 are the same Y
          if (i->getX2() < i->getX3())
          {
            v1[0] = i->getX2();
            v3[0] = i->getX3();
            v2[0] = i->getX1();
          }
          else
          {
            v1[0] = i->getX3();
            v3[0] = i->getX2();
            v2[0] = i->getX1();
          }
          if (v2[1] < v1[1])
            writeFlatTop(v1,v2,v3,screen,i);
          else
            writeFlatBottom(v1,v2,v3,screen,i);
        }
        else if (i->getY1() == i->getY3())
        {
          //1 and 3 are the same Y,
          if (i->getX1() < i->getX3())
          {
            v1[0] = i->getX1();
            v3[0] = i->getX3();
            v2[0] = i->getX2();
          }
          else
          {
            v1[0] = i->getX3();
            v3[0] = i->getX1();
            v2[0] = i->getX2();
          }
          if (v2[1] < v1[1])
            writeFlatTop(v1,v2,v3,screen,i);
          else
            writeFlatBottom(v1,v2,v3,screen,i);
        }
        else
        {
          //must split into a flat bottom and a flat top triangle

          /*Step 2*/

          i->sortVertices(v1,v2,v3);

          /*Step 3*/

          if (v1[0] == v3[0]) slope = 0; //right triangle
          else slope = (v1[1] - v3[1])/(v1[0] - v3[0]);
          intercept = v3[1] - v3[0]*slope;

          /*Step 4*/

          v4[1] = v2[1];

          if (slope == 0) v4[0] = v3[0];
          else v4[0] = (v4[1] - intercept) / slope;

          /*Step 5*/

          if (v2[0] < v4[0])
            writeFlatTop(v2,v3,v4,screen,i);
          else
            writeFlatTop(v4,v3,v2,screen,i);

          /*Step 6*/
          
          if (v2[0]<v4[0])
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
   screen.width = 1786;
   screen.height = 1344;

   // YOUR CODE GOES HERE TO DEPOSIT TRIANGLES INTO PIXELS USING THE SCANLINE ALGORITHM
   int e = 0;
  for(std::vector<Triangle>::iterator i = triangles.begin(); i != triangles.end(); ++i) 
  {
    
    depositTriangle(i,screen);
    e++;
  }
  std::cerr<<e<<" ";
  std::cerr<<"NOW IT's time To write";
   WriteImage(image, "allTriangles");
}
