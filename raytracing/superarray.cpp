/***************************************************************************
                          supergitter.cpp  -  description
                             -------------------
    begin                : Thu Dec 14 2000
    copyright            : (C) 2000 by Thomas Weigel
    email                : weigel@lat.ruhr-uni-bochum.de
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "ellipsoid.h"
#include "vector.h"
#include "superarray.h"
#include "error.h"
#include "objectshape.h"
#include "mem.h"
//#include "makros.h"
#ifdef WITH_QT
#include <qmessagebox.h>
#include <qstring.h>
#endif

int Fehler;
/*extern double Mem_usage;
extern bool use_SWAP;
extern bool canCalc;
*/

namespace GOAT
{
  namespace raytracing 
  {
    maths::Vector<std::complex<double> > INFdummy=maths::Vector<std::complex<double> > (1,1,1);

    SuperArray::SuperArray()
    {
	    H = maths::unity();
	    R = maths::unity();
	    Fehler = NO_ERRORS;
      anzEin=0;
      isequal=false;

      type=IN_HOST;
      ywerte=0;
      zwerte=0;
      G=0;
      K=0;
    }

    SuperArray::SuperArray(double r0,int nx, int ny, int nz, const int typ)
    {
	    Fehler = NO_ERRORS;
      double b=2.0*r0;
      G=0;
      isequal=false;
      anzEin=0;
      Obj=0;
      ywerte=0;
      zwerte=0;
      K=0;
      this->r0=r0;
 


      nges=maths::Vector<int>(nx,ny,nz);
      d=maths::Vector<double> (b/nx,b/ny,b/nz);
      iscleared=true;
      type=typ;
      if (type & IN_HOST)
        allockugel(); 
    }

    SuperArray::SuperArray(double r0,int nx, int ny, int nz,ObjectShape **Obj, int anzEin, const bool isAbsolute, const int typ)
    {
      Fehler = NO_ERRORS;
      double b=2.0*r0;
      isequal=isAbsolute;
      anzEin=0;
      G=0;
      this->r0=r0;
      nges=maths::Vector<int>(nx,ny,nz);
      type=typ;
      d=maths::Vector<double> (b/nx,b/ny,b/nz);
      iscleared=true;
      if(type & IN_HOST)
      {
      allockugel();
      }
      addInc(Obj,anzEin,isAbsolute);
    }


    bool SuperArray::addInc(ObjectShape **Obj,int anzEin, const bool isAbsolute)
    {
    bool ok=true;
    for (int i=0; (i<anzEin)/* && (canCalc)*/ && (ok); i++) ok=addInc (Obj[i],isAbsolute); 
    return ok;
    }

    maths::Vector<int> SuperArray::gitterpunkt(maths::Vector<double> P)
    {
    int ix,iy,iz,jx,jy,jz,jxh,jyh,anzx;
    maths::Vector<int> pi,ph;
    maths::Vector<double> h;
    h=H*P+maths::Vector<double> (r0,r0,r0);
    ph=maths::Vector<int> (floor(h[0]/d[0]),floor(h[1]/d[1]),floor(h[2]/d[2]));

    if(type & IN_HOST)
    {
      pi = ph;   
    }
    else
    {
      pi=ph;
    }
    return pi; 
    }

    bool SuperArray::addInc (ObjectShape *E, const bool isAbsolute)
    {
    int hmem,hmem2;
    char str[500];
    SysMemInfo smi;
    MemInfo mi;
    long int allocMem;
    maths::Vector<int> pul,por,N,hn;
    maths::Vector<double> h,O;
    O=maths::Vector<double> (-r0,-r0,-r0);
    
    smi=sysmem();
    mi=memstat();
      h=ceil(ediv(E->por,d))-floor(ediv(E->pul,d)); 

      hn=maths::Vector<int> ((int)h[0],(int)h[1],(int)h[2]); // Gr��e des 3D-Gitters in die drei Koordinatenrichtungen

      /* Berechne den tats�chlichen Bedarf */
      allocMem=sizeof (maths::Vector<std::complex<double> > ***)+2*sizeof (maths::Vector<int>)+sizeof (ObjectShape *)
                +hn[0]*sizeof (maths::Vector<std::complex<double> > **)
                +hn[0]*hn[1]*sizeof (maths::Vector<std::complex<double> > *)
                +hn[0]*hn[1]*hn[2]*sizeof (maths::Vector<std::complex<double> >);
    
      {
        if (anzEin<1)  // Es ist der erste Einschluss der hinzugef�gt wird
        { 
          G=(maths::Vector<std::complex<double> >  ****) malloc (sizeof (maths::Vector<std::complex<double> >  ***));
        if (G==NULL) { error (MALLOC_ERR,"SuperArray::addInc G=.."); return false; }
        Pul=(maths::Vector<int> *) malloc (sizeof (maths::Vector<int>));
      if (Pul==NULL) { error (MALLOC_ERR,"SuperArray::addInc Pul=.."); return false; }
        n=(maths::Vector<int> *) malloc (sizeof (maths::Vector<int>));
      if (n==NULL) { error (MALLOC_ERR,"SuperArray::addInc n=.."); return false; }
        Obj=(ObjectShape **) malloc (sizeof (ObjectShape *));
      if (Obj==NULL) { error (MALLOC_ERR,"SuperArray::addInc Obj=.."); return false; }
        } 
    else
      {  
        G=(maths::Vector<std::complex<double> > ****) realloc (G,(anzEin+1)*sizeof (maths::Vector<std::complex<double> >  ***));
      if (G==NULL) { error (REALLOC_ERR,"SuperArray::addInc G=.."); return false; }
        Pul=(maths::Vector<int> *) realloc (Pul,(anzEin+1)*sizeof (maths::Vector<int>));
      if (Pul==NULL) { error (REALLOC_ERR,"SuperArray::addInc G=.."); return false; }
        n=(maths::Vector<int> *) realloc (n,(anzEin+1)*sizeof (maths::Vector<int>));
      if (n==NULL) { error (REALLOC_ERR,"SuperArray::addInc n=.."); return false; }
          Obj=(ObjectShape **) realloc (Obj,(anzEin+1)*sizeof (ObjectShape *));
        if (Obj==NULL) { error (REALLOC_ERR,"SuperArray::addInc Obj=.."); return false; }
      }
      }
      
      n[anzEin]=hn;
      
      Obj[anzEin]=E;
      h=floor (ediv(Obj[anzEin]->pul+maths::Vector<double>(r0,r0,r0),d));
      Pul[anzEin]=maths::Vector<int> ((int)h[0],(int)h[1],(int)h[2]);  // ???????????
      Ellipsoid *H=(Ellipsoid *)E;
        if (E->isActive())  // Ist der Einschluss �berhaupt inelastisch aktiv ? 
        {
        G[anzEin]=(maths::Vector<std::complex<double> > ***) malloc ((n[anzEin][0]+1)*sizeof (maths::Vector<std::complex<double> > **));
        if (G[anzEin]==NULL) error(MALLOC_ERR,"SuperArray::addInc G[anzEin] G[anzEin]=..");
        for (int ix=0; ix<n[anzEin][0]+1; ix++)
        {
          G[anzEin][ix]=(maths::Vector<std::complex<double> > **) malloc ((n[anzEin][1]+1)*sizeof (maths::Vector<std::complex<double> > *));
            if (G[anzEin][ix]==NULL) error(MALLOC_ERR,"SuperArray::addInc G[anzEin][ix]=..");
            for (int iy=0; iy<n[anzEin][1]+1; iy++)
            {
              G[anzEin][ix][iy]=(maths::Vector<std::complex<double> > *) malloc ((n[anzEin][2]+1)* sizeof (maths::Vector<std::complex<double> >));
              if (G[anzEin][ix][iy]==NULL) error(MALLOC_ERR,"SuperArray::addInc G[anzEin][ix][iy]=..");
              for (int iz=0; iz<n[anzEin][2]+1; iz++)
              G[anzEin][ix][iy][iz]=maths::Vector<std::complex<double> >(0.0,0.0,0.0);
            }
        } 
        }
        else G[anzEin]=NULL;
      
      int i=0;
      anzEin++;
      iscleared=false;
      return true;
    }

    bool SuperArray::inObject (maths::Vector<double> P, int i) // da muss noch was gemacht werden
    {
    return Obj[i]->isInside (P);
    }

    bool SuperArray::inObject (maths::Vector<int> Pi, int i) 
    {
      return (Obj[i]->isInside(emult(Pi,d)));
    }

    bool SuperArray::inObject (int ix, int iy, int iz, int i)
    {
      return (Obj[i]->isInside(emult(maths::Vector<double>(ix,iy,iz),d)));
    }

    maths::Vector<std::complex<double> >& SuperArray::operator () (int ix, int iy, int iz)
    {
      maths::Vector<int> Pi=maths::Vector<int>(ix,iy,iz);
      int i=0;
      bool found=false;
      dummy=maths::Vector<std::complex<double> > (0.0,0.0,0.0); 
      if (type ==IN_OBJECT)
      {
        do
        {
          found=inObject(ix,iy,iz,i);
          i++;
        }
        while ((i<anzEin) && (!found));
        i--; 
        if (!found) return dummy;
        else
        {
          Pi=Pi-Pul[i];
          Pi=kugelindex(Pi);
          if (Fehler==NO_ERRORS)
          return G[i][Pi[0]][Pi[1]][Pi[2]]; 
          else return dummy;
        } 
      }
      else
      {
        Pi=kugelindex(Pi);
        if (Fehler==NO_ERRORS)
        return K[Pi[0]][Pi[1]][Pi[2]];
        else return dummy; 
      }
    }

    maths::Vector<std::complex<double> >& SuperArray::operator () (maths::Vector<int> Pi)
    {
      dummy=maths::Vector<std::complex<double> > (0.0,0.0,0.0); 
      if(type==IN_OBJECT)
      {
        int i=0;
        bool found=false;
        do
        {
          found=inObject(Pi,i);
          i++;
        }
        while ((i<anzEin) && (!found));
        i--;
      
        if (!found) return dummy;
        else
        {
          Pi=Pi-Pul[i];
          Pi=kugelindex(Pi);
          if (Fehler==NO_ERRORS)
          return G[i][Pi[0]][Pi[1]][Pi[2]];
          else 
          return dummy; 
        }
      }
      else
      {
        Pi=kugelindex(Pi);
        if (Fehler==NO_ERRORS)
          return K[Pi[0]][Pi[1]][Pi[2]];
        else 
          return dummy; 
      }
    }


    maths::Vector<std::complex<double> >& SuperArray::operator () (int i, int ix, int iy, int iz, bool isEinKoord)
    {
      dummy=maths::Vector<std::complex<double> >(0.0,0.0,0.0);
      maths::Vector<int> Pi=maths::Vector<int>(ix,iy,iz);
      if(type==IN_OBJECT)
      {
        if (G[i]==NULL) return dummy;
        if (!isEinKoord)  Pi=Pi-Pul[i];
        return G[i][Pi[0]][Pi[1]][Pi[2]]; 
      }
      else
      {
        Pi=kugelindex(Pi);
        if (Fehler==NO_ERRORS)
        return K[Pi[0]][Pi[1]][Pi[2]]; 
        else 
        return dummy;  
      }
    }

    maths::Vector<std::complex<double> >& SuperArray::operator () (int i, maths::Vector<int> Pi)
    {
      dummy=maths::Vector<std::complex<double> >(0.0,0.0,0.0);

      if(type==IN_OBJECT)
      {
        Pi=Pi-Pul[i];
        if (G[i]==NULL) return dummy;
        if (Pi[0]<0) return dummy; //maths::Vector<std::complex<double> > (0,0,0);
        if (Pi[1]<0) return dummy; //maths::Vector<std::complex<double> > (0,0,0);
        if (Pi[2]<0) return dummy; // maths::Vector<std::complex<double> > (0,0,0);
        if (Pi[0]>n[i][0]) 
        {
          Fehler=SUPERGITTER;
          return dummy; //maths::Vector<std::complex<double> > (0,0,0);
        }
        
        if (Pi[1]>n[i][1]) 
        {
          Fehler=SUPERGITTER;
          return dummy; //maths::Vector<std::complex<double> > (0,0,0);
        }
        
        if (Pi[2]>n[i][2]) 
        {
        Fehler=SUPERGITTER;
        }
        
        Fehler=NO_ERRORS;
        return G[i][Pi[0]][Pi[1]][Pi[2]];
      }
      else
      {
        Pi=kugelindex(Pi);
        if (Fehler==NO_ERRORS)
          return K[Pi[0]][Pi[1]][Pi[2]];
        else 
          return dummy; 
      }
    }

    maths::Vector<std::complex<double> >& SuperArray::operator () (maths::Vector<double> P)
    {
      int i;
      maths::Vector<int> Pi;
      maths::Vector<double> h;
      h=H*P-maths::Vector<double> (-r0,-r0,-r0);
      dummy=maths::Vector<std::complex<double> >(0.0,0.0,0.0);
      for (i=0; i<3; i++)
      Pi[i]=h[i]/d[i];

      if(type==IN_OBJECT)
      {
        bool found=false;
        i=0;
        do
        {
        found=inObject(Pi,i);
        i++;
        }
        while ((i<anzEin) && (!found));
        i--;
        // if (!found) return maths::Vector<std::complex<double> > (INF,INF,INF);
        if (!found) return INFdummy;
        else
        {
        if (G[i]==NULL) return dummy; 
        Pi=Pi-Pul[i];
        return G[i][Pi[0]][Pi[1]][Pi[2]];
        }
      }
      else
      {
        Pi=kugelindex(Pi);
        if (Fehler==NO_ERRORS)
        return K[Pi[0]][Pi[1]][Pi[2]];
        else return dummy;
      }
    }

    maths::Vector<std::complex<double> >& SuperArray::operator () (int i, maths::Vector<double> P)
    {
      maths::Vector<int> Pi;
      maths::Vector<double> h;
      h=H*P-maths::Vector<double> (-r0,-r0,-r0);
      dummy=maths::Vector<std::complex<double> >(0.0,0.0,0.0);

      for (int j=0; j<3; j++)
        Pi[j]=h[j]/d[j];

      if (type==IN_OBJECT)
      {
        if (G[i]==NULL) return dummy; 
        Pi=Pi-Pul[i];
        return G[i][Pi[0]][Pi[1]][Pi[2]];
      }
      else
      {
        Pi=kugelindex(Pi);
        if (Fehler==NO_ERRORS)
        return K[Pi[0]][Pi[1]][Pi[2]];
        else return dummy;
      }
    }


    void SuperArray::clear ()
    {
      int anzx, anzx2;
      anzx=nges[0];
      anzx2=nges[0]/2;

      if (type==IN_OBJECT)
      {
        if (anzEin>0)
        {
          for (int i=anzEin-1; i>=0; i--)
          {
            if (G[i]!=NULL)
            {
              for (int ix=n[i][0]; ix>=0; ix--)
              { 
                for (int iy=n[i][1]; iy>=0; iy--)
                  if (G[i][ix][iy]!=0)
                  free (G[i][ix][iy]);
                if (G[i][ix]!=0)
                {
                  free(G[i][ix]);
                } // if
              } // for ix    
            free(G[i]);
            } // if (G[i]!=NULL)
        } // for i
        free(n);
        free(Pul);
        free(Obj);
        anzEin=0;
        iscleared=true;
        G=0;
        }
      }
      else
      {
        if(K!=0)
        {
          for (int k=0;k<anzx2;k++)
          {       
            for (int l=0;l<ywerte[k];l++)
            {
              delete[] K[anzx2-1-k][ywerte[k]+l];
              delete[] K[anzx2+k][ywerte[k]+l];
              delete[] K[anzx2-1-k][ywerte[k]-1-l];
              delete[] K[anzx2+k][ywerte[k]-1-l];
            }
            delete[] K[anzx2-1-k];
            delete[] K[anzx2+k];
            delete[] zwerte[k];
          }

          delete[] ywerte;
          delete[] zwerte;

          if (anzEin>0)
          free(Obj);
        }
        K=0;
        zwerte=0;
        ywerte=0;
        anzEin=0;
        iscleared=true;
      }
    }

    void SuperArray::copy (const SuperArray &S)  // MUSS DRINGEND GE�NDERT WERDEN !
    {
    for (int i=0; i<anzEin; i++)
    if (S.G[i]!=NULL)
    {
      for (int ix=0; ix<n[i][0]; ix++)
      for (int iy=0; iy<n[i][1]; iy++)
        for (int iz=0; iz<n[i][2]; iz++)
        G[i][ix][iy][iz]=S.G[i][ix][iy][iz];
    }
}

    SuperArray& SuperArray::operator = (const SuperArray &S)
    {
      int i=0;
      int anzx, anzx2;


      if (this==&S) return *this;
      clear(); 
      r0=S.r0;
      d=S.d;
      nges=S.nges;
      type=S.type;

      anzx=nges[0]; anzx2=nges[0]/2;        

      if(type==IN_OBJECT)
      {      
        isequal=true; addInc(S.Obj,S.anzEin,true);
        isequal=true;
        if (S.anzEin!=0)
        do
        {
          for (int ix=0; ix<n[i][0]; ix++)
            for (int iy=0; iy<n[i][1]; iy++)
              for (int iz=0; iz<n[i][2]; iz++)
              G[i][ix][iy][iz]=S.G[i][ix][iy][iz];
          i++;
        } while (i<S.anzEin);
          isequal=false;
      }
      else
      {
        allockugel();
        isequal=true; addInc(S.Obj,S.anzEin,true);
        isequal=true;
        for (int k=0;k<anzx2;k++)
        {
          ywerte[k]=S.ywerte[k];
          for (int l=0;l<S.ywerte[k];l++)
          {      
            zwerte[k][l]=S.zwerte[k][l];
            for (int m=0;m<S.zwerte[k][l];m++)
            {
              K[anzx2-1-k][ywerte[k]-1-l][m]=S.K[anzx2-1-k][ywerte[k]-1-l][m];
              K[anzx2+k][ywerte[k]-1-l][m] = S.K[anzx2+k][ywerte[k]-1-l][m];
              K[anzx2-1-k][ywerte[k]+l][m]= S.K[anzx2-1-k][ywerte[k]+l][m];
              K[anzx2+k][ywerte[k]+l][m] = S.K[anzx2+k][ywerte[k]+l][m];
              K[anzx2-1-k][ywerte[k]-1-l][2*zwerte[k][l]-1-m]= S.K[anzx2-1-k][ywerte[k]-1-l][2*zwerte[k][l]-1-m];
              K[anzx2+k][ywerte[k]-1-l][2*zwerte[k][l]-1-m] = S.K[anzx2+k][ywerte[k]-1-l][2*zwerte[k][l]-1-m];
              K[anzx2-1-k][ywerte[k]+l][2*zwerte[k][l]-1-m]= S.K[anzx2-1-k][ywerte[k]+l][2*zwerte[k][l]-1-m];
              K[anzx2+k][ywerte[k]+l][2*zwerte[k][l]-1-m] = S.K[anzx2+k][ywerte[k]+l][2*zwerte[k][l]-1-m];
            }
          }
        }
      }
      return *this;
    }



    void SuperArray::add (const SuperArray& S)
    {
      int anzx2=nges[0]/2;
      isequal=false;
      H = S.H;
      R = S.R;
      if(type==IN_OBJECT)
      {
      for (int i=0; i<S.anzEin; i++)
        for (int ix=0; ix<S.n[i][0]; ix++)
          for (int iy=0; iy<S.n[i][1]; iy++)
            for (int iz=0; iz<S.n[i][2]; iz++)
            G[i][ix][iy][iz]+=S.G[i][ix][iy][iz];
      }
      else
      {
        for (int k=0;k<anzx2;k++)
        {
          for (int l=0;l<ywerte[k];l++)
          {

            for (int m=0;m<zwerte[k][l];m++)
            {
              K[anzx2-1-k][ywerte[k]-1-l][m]+=S.K[anzx2-1-k][ywerte[k]-1-l][m];
              K[anzx2+k][ywerte[k]-1-l][m] += S.K[anzx2+k][ywerte[k]-1-l][m];
              K[anzx2-1-k][ywerte[k]+l][m]+=S.K[anzx2-1-k][ywerte[k]+l][m];
              K[anzx2+k][ywerte[k]+l][m] +=S.K[anzx2+k][ywerte[k]+l][m];
              K[anzx2-1-k][ywerte[k]-1-l][2*zwerte[k][l]-1-m]+=S.K[anzx2-1-k][ywerte[k]-1-l][2*zwerte[k][l]-1-m];
              K[anzx2+k][ywerte[k]-1-l][2*zwerte[k][l]-1-m] +=S.K[anzx2+k][ywerte[k]-1-l][2*zwerte[k][l]-1-m];
              K[anzx2-1-k][ywerte[k]+l][2*zwerte[k][l]-1-m]+=S.K[anzx2-1-k][ywerte[k]+l][2*zwerte[k][l]-1-m];
              K[anzx2+k][ywerte[k]+l][2*zwerte[k][l]-1-m] +=S.K[anzx2+k][ywerte[k]+l][2*zwerte[k][l]-1-m];
            }
          }
        }
      }
    }

    void SuperArray::sub (const SuperArray& S)
    {
      int anzx2=nges[0]/2;
      
      isequal=false;
      if(type==IN_OBJECT)
      {
        for (int i=0; i<S.anzEin; i++)
          for (int ix=0; ix<S.n[i][0]; ix++)
            for (int iy=0; iy<S.n[i][1]; iy++)
              for (int iz=0; iz<S.n[i][2]; iz++)
              G[i][ix][iy][iz]-=S.G[i][ix][iy][iz];
      }
      else
      {
        for (int k=0;k<anzx2;k++)
        {
          for (int l=0;l<ywerte[k];l++)
          { 

            for (int m=0;m<zwerte[k][l];m++)
            {
              K[anzx2-1-k][ywerte[k]-1-l][m]-=S.K[anzx2-1-k][ywerte[k]-1-l][m];
              K[anzx2+k][ywerte[k]-1-l][m] -= S.K[anzx2+k][ywerte[k]-1-l][m];
              K[anzx2-1-k][ywerte[k]+l][m]-=S.K[anzx2-1-k][ywerte[k]+l][m];
              K[anzx2+k][ywerte[k]+l][m] -=S.K[anzx2+k][ywerte[k]+l][m];
              K[anzx2-1-k][ywerte[k]-1-l][2*zwerte[k][l]-1-m]-=S.K[anzx2-1-k][ywerte[k]-1-l][2*zwerte[k][l]-1-m];
              K[anzx2+k][ywerte[k]-1-l][2*zwerte[k][l]-1-m] -=S.K[anzx2+k][ywerte[k]-1-l][2*zwerte[k][l]-1-m];
              K[anzx2-1-k][ywerte[k]+l][2*zwerte[k][l]-1-m]-=S.K[anzx2-1-k][ywerte[k]+l][2*zwerte[k][l]-1-m];
              K[anzx2+k][ywerte[k]+l][2*zwerte[k][l]-1-m] -=S.K[anzx2+k][ywerte[k]+l][2*zwerte[k][l]-1-m];
            }
          }
        }
      }
    }

    std::ostream&   operator << (std::ostream &os, const SuperArray &S)
    {
      os << "r0=" << S.r0 << std::endl;
      os << "Ausdehnung:    nx=" << S.nges[0] << "  ny=" << S.nges[1] << "  nz=" << S.nges[2] << std::endl;
      os << S.anzEin << " Einschluesse" << std::endl;
      if (S.anzEin>0)
      for (int i=0; i<S.anzEin; i++)
      {
        os << "====================== Einschluss Nr. " << i << " ======================" << std::endl;
        for (int ix=0; ix<S.n[i][0]; ix++)
          for (int iy=0; iy<S.n[i][1]; iy++)
            for (int iz=0; iz<S.n[i][2]; iz++)
              os <<  abs(S.G[i][ix][iy][iz]) << std::endl;
      }
      return os;
    }

    void SuperArray::fill (const maths::Vector<std::complex<double> > &x)
    {
      int anzx,anzx2;
      anzx=nges[0]; anzx2=nges[0]/2;

      if(type==IN_OBJECT)
      {
      for  (int i=0; i<anzEin; i++)
        for (int ix=0; ix<n[i][0]; ix++)
        for (int iy=0; iy<n[i][1]; iy++)
          for (int iz=0; iz<n[i][2]; iz++)
            G[i][ix][iy][iz]=x;
      }
      else
      {
        for (int k=0;k<anzx2;k++)
        { 
          for (int l=0;l<ywerte[k];l++)
          {
            for (int m=0;m<zwerte[k][l];m++)
            {
              K[anzx2-1-k][ywerte[k]-1-l][zwerte[k][l]-1-m]=x;
              K[anzx2+k][ywerte[k]-1-l][zwerte[k][l]-1-m] = x;
              K[anzx2-1-k][ywerte[k]+l][zwerte[k][l]-1-m]=x;
              K[anzx2+k][ywerte[k]+l][zwerte[k][l]-1-m] =x;
              K[anzx2-1-k][ywerte[k]-1-l][zwerte[k][l]+m]=x;
              K[anzx2+k][ywerte[k]-1-l][zwerte[k][l]+m] =x;
              K[anzx2-1-k][ywerte[k]+l][zwerte[k][l]+m]=x;
              K[anzx2+k][ywerte[k]+l][zwerte[k][l]+m] =x;
            }
          }
        }
      }
    }

 
    void SuperArray::makeReal ()
    {
    for (int i=0; i<anzEin; i++)
      for (int ix=0; ix<n[i][0]; ix++)
       for (int iy=0; iy<n[i][1]; iy++)
        for (int iz=0; iz<n[i][2]; iz++)
          for (int j=0; j<3; j++)
          G[i][ix][iy][iz][j]=abs(G[i][ix][iy][iz][j]);
    }



    void SuperArray::saveExPhase (char* FName,int i)
    {
      maths::Vector<int> Pi;
      maths::Vector<std::complex<double> > E;
      std::ofstream os;
      os.open (FName);
      std::complex<double> phase;
      if (type==IN_OBJECT)
      {
        os<<"%Dimensionen "<<n[i][0]<<"  x  "<<  n[i][1]<<"  x  "<<n[i][2]<< std::endl;

        for (int ix=0; ix<n[i][0]; ix++)
        for (int iy=0; iy<n[i][1]; iy++)
          for (int iz=0; iz<n[i][2]; iz++)
          {       
            os << atan2(imag(G[i][ix][iy][iz][0]),real(G[i][ix][iy][iz][0])) << std::endl;     
          }   
      }
      else
      {
        for (int ix=0; ix<nges[0];ix++)
          for (int iy=0; iy<nges[1];iy++)
          for (int iz=0; iz<nges[2];iz++)
          {                     
            Pi=kugelindex(maths::Vector<int>(ix,iy,iz));
            if (Fehler!=NO_ERRORS) os << 0.0 << std::endl;
            else
            {
              E=K[Pi[0]][Pi[1]][Pi[2]];	
              os << atan2(imag(E[0]),real(E[0])) << std::endl;	
            }
          }
      }
      os.close();
    }

void SuperArray::saveEyPhase (char* FName,int i)
{
  maths::Vector<int> Pi;
  maths::Vector<std::complex<double> > E;
  std::ofstream os;
  os.open (FName);
  std::complex<double> phase;
if (type==IN_OBJECT)
{
  os<<"%Dimensionen "<<n[i][0]<<"  x  "<<  n[i][1]<<"  x  "<<n[i][2]<< std::endl;

  for (int ix=0; ix<n[i][0]; ix++)
   for (int iy=0; iy<n[i][1]; iy++)
    for (int iz=0; iz<n[i][2]; iz++)
    {
       os << atan2(imag(G[i][ix][iy][iz][1]),real(G[i][ix][iy][iz][1])) << std::endl;             
    }   
}
else
{
   for (int ix=0; ix<nges[0];ix++)
    for (int iy=0; iy<nges[1];iy++)
     for (int iz=0; iz<nges[2];iz++)
     {                     
       Pi=kugelindex(maths::Vector<int>(ix,iy,iz));
       if (Fehler!=NO_ERRORS) os << 0.0 << std::endl;
       else
       {
        E=K[Pi[0]][Pi[1]][Pi[2]];
	os << atan2(imag(E[1]),real(E[1])) << std::endl;	
       }
     }
}
  os.close();
}

void SuperArray::saveEzPhase (char* FName,int i)
{
  maths::Vector<int> Pi;
  maths::Vector<std::complex<double> > E;
  std::ofstream os;
  os.open (FName);
  std::complex<double> phase;
if (type==IN_OBJECT)
{
  os<<"%Dimensionen "<<n[i][0]<<"  x  "<<  n[i][1]<<"  x  "<<n[i][2]<<std::endl;

  for (int ix=0; ix<n[i][0]; ix++)
   for (int iy=0; iy<n[i][1]; iy++)
    for (int iz=0; iz<n[i][2]; iz++)
    {       
       os << atan2(imag(G[i][ix][iy][iz][2]),real(G[i][ix][iy][iz][2])) << std::endl;
    }   
}
else
{
   for (int ix=0; ix<nges[0];ix++)
    for (int iy=0; iy<nges[1];iy++)
     for (int iz=0; iz<nges[2];iz++)
     {                     
       Pi=kugelindex(maths::Vector<int>(ix,iy,iz));
       if (Fehler!=NO_ERRORS) os << 0.0 << std::endl;
       else
       {
        E=K[Pi[0]][Pi[1]][Pi[2]];
	os << atan2(imag(E[2]),real(E[2])) << std::endl;
       }
     }
}
  os.close();
}


void SuperArray::saveExPol (char* FName,int i)
{
  maths::Vector<int> Pi;
  maths::Vector<std::complex<double> > E;
  std::ofstream os;
  os.open (FName);
if (type==IN_OBJECT)
{
  os<<"%Dimensionen "<<n[i][0]<<"  x  "<<  n[i][1]<<"  x  "<<n[i][2]<<std::endl;

  for (int ix=0; ix<n[i][0]; ix++)
   for (int iy=0; iy<n[i][1]; iy++)
    for (int iz=0; iz<n[i][2]; iz++)
       os << abs(G[i][ix][iy][iz][0]) << std::endl;
}
else
{
   for (int ix=0; ix<nges[0];ix++)
    for (int iy=0; iy<nges[1];iy++)
     for (int iz=0; iz<nges[2];iz++)
     {                     
       Pi=kugelindex(maths::Vector<int>(ix,iy,iz));
       if (Fehler!=NO_ERRORS) os << 0.0 << std::endl;
       else
       {
        E=K[Pi[0]][Pi[1]][Pi[2]];
        os << abs(E[0]) << std::endl;
       }
     }
}
  os.close();
}

void SuperArray::saveEyPol (char* FName,int i)
{
  maths::Vector<int> Pi;
  maths::Vector<std::complex<double> > E;
  std::ofstream os;
  os.open (FName);
if (type==IN_OBJECT)
{
  os<<"%Dimension "<<n[i][0]<<"  x  "<<  n[i][1]<<"  x  "<<n[i][2]<<std::endl;
  for (int ix=0; ix<n[i][0]; ix++)
   for (int iy=0; iy<n[i][1]; iy++)
    for (int iz=0; iz<n[i][2]; iz++)
     // os << abs(G[i][ix][iy][iz][1])/abs(G[i][ix][iy][iz]) << std::endl;
     //  os << real(G[i][ix][iy][iz][1]) << std::endl;
      os << abs(G[i][ix][iy][iz][1]) << std::endl;
  }
  else
  {
   for (int ix=0; ix<nges[0];ix++)
    for (int iy=0; iy<nges[1];iy++)
     for (int iz=0; iz<nges[2];iz++)
     {
       Pi=kugelindex(maths::Vector<int>(ix,iy,iz));
       if (Fehler!=NO_ERRORS) os << 0.0 << std::endl;
       else
       {
        E=K[Pi[0]][Pi[1]][Pi[2]];
        os << abs(E[1]) << std::endl;
       }
     }


  }

  
  os.close();
}

void SuperArray::saveEzPol (char* FName,int i)
{
  maths::Vector<int> Pi;
  maths::Vector<std::complex<double> > E;
  std::ofstream os;
  os.open (FName);
if (type==IN_OBJECT)
{
  os<<"%Dimension "<<n[i][0]<<"  x  "<<  n[i][1]<<"  x  "<<n[i][2]<<std::endl;


  for (int ix=0; ix<n[i][0]; ix++)
   for (int iy=0; iy<n[i][1]; iy++)
    for (int iz=0; iz<n[i][2]; iz++)
  //   os << abs(G[i][ix][iy][iz][2])/abs(G[i][ix][iy][iz]) << std::endl;
  // os << real(G[i][ix][iy][iz][2]) << std::endl;
  os << abs(G[i][ix][iy][iz][2]) << std::endl;
}
else
{
   for (int ix=0; ix<nges[0];ix++)
    for (int iy=0; iy<nges[1];iy++)
     for (int iz=0; iz<nges[2];iz++)
     {
       Pi=kugelindex(maths::Vector<int>(ix,iy,iz));
       if (Fehler!=NO_ERRORS) os << 0.0 << std::endl;
       else
       {
        E=K[Pi[0]][Pi[1]][Pi[2]];
        os << abs(E[2]) << std::endl;
       }
     }

}

  os.close();
}

void SuperArray::saveabsE (const char* FName,int i)
{
  maths::Vector<std::complex<double> > E;
  std::ofstream os;
  os.open (FName);
  maths::Vector<int> Pi;
  double x;
  if(type==IN_OBJECT)
  {
  os<<"%Dimension "<<n[i][0]<<"  x  "<<  n[i][1]<<"  x  "<<n[i][2]<<std::endl;


  for (int ix=0; ix<n[i][0]; ix++)
   for (int iy=0; iy<n[i][1]; iy++)
    for (int iz=0; iz<n[i][2]; iz++)
    {
      x=real(G[i][ix][iy][iz][0]*conj(G[i][ix][iy][iz][0]));
      x+=real(G[i][ix][iy][iz][1]*conj(G[i][ix][iy][iz][1]));
      x+=real(G[i][ix][iy][iz][2]*conj(G[i][ix][iy][iz][2]));
      os << x << std::endl;
    }
  }
  else
 {
   for (int ix=0; ix<nges[0];ix++)
    for (int iy=0; iy<nges[1];iy++)
     for (int iz=0; iz<nges[2];iz++)
     {
       Pi=kugelindex(maths::Vector<int>(ix,iy,iz));
       if (Fehler!=NO_ERRORS) os << 0.0 << std::endl;
       else 
       {
        E=K[Pi[0]][Pi[1]][Pi[2]];
        os << real(E*conj(E)) << std::endl; 
       } 
     }
 }
  
  os.close();
}

void SuperArray::saveFullE(const char* FName, int i)
{
    maths::Vector<std::complex<double> > E;
    std::ofstream os;
    os.open(FName);
    maths::Vector<int> Pi;
    double x;
    if (type == IN_OBJECT)
    {
        os << "%Dimension " << n[i][0] << "  x  " << n[i][1] << "  x  " << n[i][2] << std::endl;


        for (int ix = 0; ix < n[i][0]; ix++)
            for (int iy = 0; iy < n[i][1]; iy++)
                for (int iz = 0; iz < n[i][2]; iz++)
                {
                    os << real(G[i][ix][iy][iz][0]) << "\t" << imag(G[i][ix][iy][iz][0]) << "\t";
                    os << real(G[i][ix][iy][iz][1]) << "\t" << imag(G[i][ix][iy][iz][1]) << "\t";
                    os << real(G[i][ix][iy][iz][2]) << "\t" << imag(G[i][ix][iy][iz][2]) << "\t" << std::endl;
                }
    }
    else
    {
        for (int ix = 0; ix < nges[0]; ix++)
            for (int iy = 0; iy < nges[1]; iy++)
                for (int iz = 0; iz < nges[2]; iz++)
                {
                    Pi = kugelindex(maths::Vector<int>(ix, iy, iz));
                    if (Fehler != NO_ERRORS) os << 0.0 << std::endl;
                    else
                    {
                        E = K[Pi[0]][Pi[1]][Pi[2]];
                        os << real(E[0]) << "\t" << imag(E[0]) << "\t";
                        os << real(E[1]) << "\t" << imag(E[1]) << "\t";
                        os << real(E[2]) << "\t" << imag(E[2]) << "\t" << std::endl;
                    }
                }
    }

    os.close();
}




double sumabs (const SuperArray &S)
{
  double Erg,h;
  int anzx2=S.nges[0]/2;
  Erg=0.0;
 if(S.type==IN_OBJECT)
 {
  for (int i=0; i<S.anzEin; i++)
  for (int ix=0; ix<S.n[i][0]; ix++)
   for (int iy=0; iy<S.n[i][1]; iy++)
    for (int iz=0; iz<S.n[i][2]; iz++)
     {
       h=abs(S.G[i][ix][iy][iz]);
      Erg+=h;
     }
 }
 else
 {
  for (int k=0;k<anzx2;k++)
  {
    for (int l=0;l<S.ywerte[k];l++)
    {
 
      for (int m=0;m<S.zwerte[k][l];m++)
      {
        h=abs(S.K[anzx2-1-k][S.ywerte[k]-1-l][m]);
        Erg+=h;
        h=abs(S.K[anzx2+k][S.ywerte[k]-1-l][m]);
        Erg+=h;
        h=abs(S.K[anzx2-1-k][S.ywerte[k]+l][m]);
        Erg+=h;
        h=abs(S.K[anzx2+k][S.ywerte[k]+l][m]);
        Erg+=h;
        h=abs(S.K[anzx2-1-k][S.ywerte[k]-1-l][2*S.zwerte[k][l]-1-m]);
        Erg+=h;
        h=abs(S.K[anzx2+k][S.ywerte[k]-1-l][2*S.zwerte[k][l]-1-m]);
        Erg+=h;
        h=abs(S.K[anzx2-1-k][S.ywerte[k]+l][2*S.zwerte[k][l]-1-m]);
        Erg+=h;
        h=abs(S.K[anzx2+k][S.ywerte[k]+l][2*S.zwerte[k][l]-1-m]);
        Erg+=h;
      }
    }
  }
  
 }
 
  return Erg;
}

double sumabs2 (const SuperArray &S)
{
  double Erg,h;
  int anzx2=S.nges[0]/2;
  Erg=0.0;
 if(S.type==IN_OBJECT)
 {
  for (int i=0; i<S.anzEin; i++)
  for (int ix=0; ix<S.n[i][0]; ix++)
   for (int iy=0; iy<S.n[i][1]; iy++)
    for (int iz=0; iz<S.n[i][2]; iz++)
     {
       h=abs(S.G[i][ix][iy][iz]);
       Erg+=h*h;
     }
 }
 else
 {
  for (int k=0;k<anzx2;k++)
  {
    for (int l=0;l<S.ywerte[k];l++)
    {
 
      for (int m=0;m<S.zwerte[k][l];m++)
      {
        h=abs(S.K[anzx2-1-k][S.ywerte[k]-1-l][m]);
        Erg+=h*h;
        h=abs(S.K[anzx2+k][S.ywerte[k]-1-l][m]);
        Erg+=h*h;
        h=abs(S.K[anzx2-1-k][S.ywerte[k]+l][m]);
        Erg+=h*h;
        h=abs(S.K[anzx2+k][S.ywerte[k]+l][m]);
        Erg+=h*h;
        h=abs(S.K[anzx2-1-k][S.ywerte[k]-1-l][2*S.zwerte[k][l]-1-m]);
        Erg+=h*h;
        h=abs(S.K[anzx2+k][S.ywerte[k]-1-l][2*S.zwerte[k][l]-1-m]);
        Erg+=h*h;
        h=abs(S.K[anzx2-1-k][S.ywerte[k]+l][2*S.zwerte[k][l]-1-m]);
        Erg+=h*h;
        h=abs(S.K[anzx2+k][S.ywerte[k]+l][2*S.zwerte[k][l]-1-m]);
        Erg+=h*h;
      }
    }
  }

 }
  return Erg;
}

double abs2sum (const SuperArray &S)
{
 maths::Vector<std::complex<double> > h;
 double Erg;
 int anzx2=S.nges[0]/2;
 if(S.type==IN_OBJECT)
 {
 for (int i=0; i<S.anzEin; i++)
  for (int ix=0; ix<S.n[i][0]; ix++)
   for (int iy=0; iy<S.n[i][1]; iy++)
    for (int iz=0; iz<S.n[i][2]; iz++)
     {
       h+=S.G[i][ix][iy][iz];
     }
 }
 else
 {
    for (int k=0;k<anzx2;k++)
  {
    for (int l=0;l<S.ywerte[k];l++)
    {

      for (int m=0;m<S.zwerte[k][l];m++)
      {
        h+=S.K[anzx2-1-k][S.ywerte[k]-1-l][m];
        h+=S.K[anzx2+k][S.ywerte[k]-1-l][m];
        h+=S.K[anzx2-1-k][S.ywerte[k]+l][m];
        h+=S.K[anzx2+k][S.ywerte[k]+l][m];
        h+=S.K[anzx2-1-k][S.ywerte[k]-1-l][2*S.zwerte[k][l]-1-m];
        h+=S.K[anzx2+k][S.ywerte[k]-1-l][2*S.zwerte[k][l]-1-m];
        h+=S.K[anzx2-1-k][S.ywerte[k]+l][2*S.zwerte[k][l]-1-m];
        h+=S.K[anzx2+k][S.ywerte[k]+l][2*S.zwerte[k][l]-1-m];
      }
    }
  }

 }
 Erg=abs(h);
 return Erg*Erg;
}

void SuperArray::allockugel ()
{
  int anzx, anzx2; 
  double dx,dy,dz,hilf;
  anzx=nges[0];
  anzx2=nges[0]/2;
 // cout << "allockugel" << std::endl;
  ywerte = new int[anzx2];
  zwerte = new int*[anzx2];

  /*dx=d[0];
  dy=d[1];
  dz=d[2];*/
  dx=2.0/nges[0];
  dy=2.0/nges[1];
  dz=2.0/nges[2];

  // cout << "dx=" << dx << std::endl; 
  for(int i=0;i<anzx2;i++)
  {
    ywerte[i]=int(ceil(sqrt(1.0-(i*dx)*(i*dx))/dy));
    zwerte[i]=new int[anzx2];
    for(int j=0;j<anzx2;j++)
    {
      hilf=1.0-(i*dx)*(i*dx)-(j*dy)*(j*dy);
//      cout << "hilf:" << hilf << std::endl;
      if (hilf>=0)
//        zwerte[i][j]=int(ceil(sqrt(1.0-(i*dx)*(i*dx)-(j*dy)*(j*dy))*nges[2]/2.0));
          zwerte[i][j]=ceil(sqrt(hilf)/dz);
      else
        zwerte[i][j]=1;

//      cout << "zwerte[" << i << "][" << j << "]:" <<  zwerte[i][j] << std::endl;
    }
  }

 K=new maths::Vector<std::complex<double> >**[anzx];
 for (int k=0;k<anzx2;k++)
 {
   K[anzx2-1-k]=new maths::Vector<std::complex<double> >*[2*ywerte[k]];
   K[anzx2+k]=new maths::Vector<std::complex<double> >*[2*ywerte[k]];
   for (int l=0;l<ywerte[k];l++)
   {
     K[anzx2-1-k][ywerte[k]+l]=new maths::Vector<std::complex<double> >[2*zwerte[k][l]];
//    cout << "K[" << anzx2-1-k << "][" << ywerte[k]+l<<"]=" << 2*zwerte[k][l] << std::endl;

     K[anzx2+k][ywerte[k]+l]=new maths::Vector<std::complex<double> >[2*zwerte[k][l]];
//    cout << "K[" << anzx2+k << "][" << ywerte[k]+l<<"]=" << 2*zwerte[k][l] << std::endl;

     K[anzx2-1-k][ywerte[k]-1-l]=new maths::Vector<std::complex<double> >[2*zwerte[k][l]];
//    cout << "K[" << anzx2-1-k << "][" << ywerte[k]-1-l<<"]=" << 2*zwerte[k][l] << std::endl;

     K[anzx2+k][ywerte[k]-1-l]=new maths::Vector<std::complex<double> >[2*zwerte[k][l]];
//     cout << "K[" << anzx2+k << "][" << ywerte[k]-1-l<<"]=" << zwerte[k][l] << std::endl;
   }
 }
}

maths::Vector<int> SuperArray::kugelindex(maths::Vector<int> Pi)
{
 int anzx,jx,jy,jz,ix,iy,iz,jxh,jyh;
 maths::Vector<int> pk;
 anzx=nges[0];
 ix=Pi[0]; iy=Pi[1]; iz=Pi[2]; 
 Fehler=NO_ERRORS;
 jx=ix;
 jxh=abs(int(jx-(anzx-1)/2.0));
 jy=iy-(anzx/2-ywerte[jxh]);
 if ((jy<0)||(jy>(2*ywerte[jxh]-1)))   // Punkt ausserhalb 
 { 
    Fehler=SUPERGITTER;
//    cout << "ix:" << ix << ", iy:" << iy << ", iz:" << iz <<  ", jx:" << jx <<
//             ", jxh:" << jxh << ", jy:" << jy << std::endl;
//    cout << "ywerte[" << jxh << "]:" << ywerte[jxh] << std::endl;
    return maths::Vector<int> (-1,-1,-1);
 } 
 else
 {
  jyh=abs(int(jy-(2*ywerte[jxh]-1)/2.0));
  jz=iz-(anzx/2-zwerte[jxh][jyh]);

  if ((jz<0)||(jz>(2*zwerte[jxh][jyh])-1)) // Punkt ausserhalb 
  {
   Fehler=SUPERGITTER;
//    cout << "ix:" << ix << ", iy:" << iy << ", iz:" << iz << ", jx:" << jx << ", jy:" << jy << 
//             ", jxh:" << jxh << ", jyh:" << jyh << ", jz:" << jz << std::endl;
//    cout << "ywerte[" << jxh << "]:" << ywerte[jxh] << ", zwerte[" << jxh << "][" << jyh << "]:" <<
//    zwerte[jxh][jyh] << std::endl;
    return maths::Vector<int> (-1,-1,-1);
   return maths::Vector<int> (-1,-1,-1);
  }
  else
  {
   return  maths::Vector<int>(jx,jy,jz);
  }
 }
 return pk;
}

maths::Vector<std::complex<double> > SuperArray::kugelwert(maths::Vector<int> Pi)
{
 int anzx,jx,jy,jz,ix,iy,iz,jxh,jyh;
// maths::Vector<std::complex<double> > pc;
 anzx=nges[0];
 ix=Pi[0]; iy=Pi[1]; iz=Pi[2];
 jx=ix;
 jxh=abs(int(jx-(anzx-1)/2.0));
 jy=iy-(anzx/2-ywerte[jxh]);
 pc=maths::Vector<std::complex<double> >(0.0,0.0,0.0);
 if ((jy<0)||(jy>(2*ywerte[jxh]-1)))
 {
  std::cout << "Fehler!! Index iy falsch" << std::endl;
    return pc;
 }
 else
 {
  jyh=abs(int(jy-(2*ywerte[jxh]-1)/2.0));
  jz=iz-(anzx/2-zwerte[jxh][jyh]);

  if ((jz<0)||(jz>(2*zwerte[jxh][jyh])-1))
  {
   return pc;
  }
  return K[jx][jy][jz];
 }
 return pc;
}

maths::Vector<std::complex<double> > SuperArray::kugelwert(int ix, int iy, int iz)
{
 int anzx,jx,jy,jz,jxh,jyh;
 maths::Vector<std::complex<double> > pc;
 anzx=nges[0];


 jx=ix;
 jxh=abs(int(jx-(anzx-1)/2.0));
 jy=iy-(anzx/2-ywerte[jxh]);
 if ((jy<0)||(jy>(2*ywerte[jxh]-1)))
 {
  std::cout << "Fehler!! Index iy falsch" << std::endl;
  pc=maths::Vector<std::complex<double> >(0.0,0.0,0.0);
  return pc;
 }
 else
 {
  jyh=abs(int(jy-(2*ywerte[jxh]-1)/2.0));
  jz=iz-(anzx/2-zwerte[jxh][jyh]);
  if ((jz<0)||(jz>(2*zwerte[jxh][jyh])-1))
  {
   std::cout << "Fehler!! Index iz falsch" << std::endl;
   pc=maths::Vector<std::complex<double> >(0.0,0.0,0.0);
   return pc;
  }
  pc=K[jx][jy][jz];
 }
 return pc;
}
  }
}