#include "matrix.h"

using namespace std;
namespace GOAT
{
    namespace maths
    {
        Matrix<complex<double> > operator * (const Matrix<double>& A, const Matrix<complex<double> >& B)
        {
            Matrix<complex<double> > Erg;
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    for (int l = 0; l < 3; l++)
                        Erg.M[i][j] += A.M[i][l] * B.M[l][j];
            return Erg;
        }

        Matrix<complex<double> > operator * (const Matrix<complex<double> >& A, const Matrix<double>& B)
        {
            Matrix<complex<double> > Erg;
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    for (int l = 0; l < 3; l++)
                        Erg.M[i][j] += A.M[i][l] * B.M[l][j];
            return Erg;
        }

        Matrix<complex<double> > operator + (const Matrix<double>& A, const Matrix<complex<double> >& B)
        {
            Matrix<complex<double> > Erg;
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    Erg.M[i][j] = A.M[i][j] + B.M[i][j];
            return Erg;
        }

        Matrix<complex<double> > operator + (const Matrix<complex<double> >& A, const Matrix<double>& B)
        {
            Matrix<complex<double> > Erg;
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    Erg.M[i][j] = A.M[i][j] + B.M[i][j];
            return Erg;
        }

        Matrix<complex<double> > operator - (const Matrix<double>& A, const Matrix<complex<double> >& B)
        {
            Matrix<complex<double> > Erg;
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    Erg.M[i][j] = A.M[i][j] - B.M[i][j];
            return Erg;
        }

        Matrix<complex<double> > operator - (const Matrix<complex<double> >& A, const Matrix<double>& B)
        {
            Matrix<complex<double> > Erg;
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    Erg.M[i][j] = A.M[i][j] - B.M[i][j];
            return Erg;
        }

        Vector<complex<double> > operator * (const Matrix<double>& A, const Vector<complex<double> >& r)
        {
            Vector<complex<double> > Erg;
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    Erg[i] += A.M[i][j] * r[j];
            return Erg;
        }

        Vector<complex<double> > operator * (const Matrix<complex<double> >& A, const Vector<double>& r)
        {
            Vector<complex<double> > Erg;
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    Erg[i] += A.M[i][j] * r[j];
            return Erg;
        }

        Matrix<complex<double> > operator * (const Matrix<double>& A, const complex<double>& x)
        {
            Matrix<complex<double> > Erg;
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    Erg.M[i][j] = A.M[i][j] * x;
            return Erg;
        }

        Matrix<complex<double> > operator * (const Matrix<complex<double> >& A, const double& x)
        {
            Matrix<complex<double> > Erg;
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    Erg.M[i][j] = A.M[i][j] * x;
            return Erg;
        }

        Matrix<complex<double> > operator * (const double& x, const Matrix<complex<double> >& A)
        {
            Matrix<complex<double> > Erg;
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    Erg.M[i][j] = A.M[i][j] * x;
            return Erg;
        }

        Matrix<complex<double> > operator * (const complex<double>& x, const Matrix<double>& A)
        {
            Matrix<complex<double> > Erg;
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    Erg.M[i][j] = A.M[i][j] * x;
            return Erg;
        }

        Matrix<complex<double> > operator / (const Matrix<double>& A, const complex<double>& x)
        {
            Matrix<complex<double> > Erg;
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    Erg.M[i][j] = A.M[i][j] / x;
            return Erg;
        }

        Matrix<complex<double> > operator / (const Matrix<complex<double> >& A, const double& x)
        {
            Matrix<complex<double> > Erg;
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    Erg.M[i][j] = A.M[i][j] / x;
            return Erg;
        }



        Matrix<double> unity()
        {
            Matrix<double> Erg;
            Erg.M[0][0] = 1.0;  Erg.M[0][1] = 0.0; Erg.M[0][2] = 0.0;
            Erg.M[1][0] = 0.0;  Erg.M[1][1] = 1.0; Erg.M[1][2] = 0.0;
            Erg.M[2][0] = 0.0;  Erg.M[2][1] = 0.0; Erg.M[2][2] = 1.0;
            return Erg;
        }
        Matrix<complex<double> > cunity()
        {
            Matrix<complex<double> > Erg;
            Erg.M[0][0] = 1.0;  Erg.M[0][1] = 0.0; Erg.M[0][2] = 0.0;
            Erg.M[1][0] = 0.0;  Erg.M[1][1] = 1.0; Erg.M[1][2] = 0.0;
            Erg.M[2][0] = 0.0;  Erg.M[2][1] = 0.0; Erg.M[2][2] = 1.0;
            return Erg;
        }


        Matrix<double> null()
        {
            Matrix<double> Erg;
            Erg.M[0][0] = 0.0;  Erg.M[0][1] = 0.0; Erg.M[0][2] = 0.0;
            Erg.M[1][0] = 0.0;  Erg.M[1][1] = 0.0; Erg.M[1][2] = 0.0;
            Erg.M[2][0] = 0.0;  Erg.M[2][1] = 0.0; Erg.M[2][2] = 0.0;
            return Erg;
        }

        Matrix<double> drehz(const double& gamma)
        {
            Matrix<double> M;
            M.M[0][0] = cos(gamma); M.M[0][1] = -sin(gamma); M.M[0][2] = 0.0;
            M.M[1][0] = sin(gamma); M.M[1][1] = cos(gamma); M.M[1][2] = 0.0;
            M.M[2][0] = 0.0; M.M[2][1] = 0.0; M.M[2][2] = 1.0;
            return M;
        }

        Matrix<double> rotMatrix(const Vector<double> a, double gamma)
            /*
              erstellt die rotMatrix für eine Drehung um einen Winkel gamma
              die Drehachse wird durch den Vektor a dargestellt
              (Drehrichtung: gegen den Uhrzeigersinn)
            */

        {
            Matrix<double> A, B, Erg;
            A(0, 0) = a[0] * a[0];    A(0, 1) = a[0] * a[1];    A(0, 2) = a[0] * a[2];
            A(1, 0) = a[1] * a[0];    A(1, 1) = a[1] * a[1];    A(1, 2) = a[1] * a[2];
            A(2, 0) = a[2] * a[0];    A(2, 1) = a[2] * a[1];    A(2, 2) = a[2] * a[2];

            B(0, 0) = 0.0;    B(0, 1) = -a[2];    B(0, 2) = a[1];
            B(1, 0) = a[2];   B(1, 1) = 0.0;      B(1, 2) = -a[0];
            B(2, 0) = -a[1];  B(2, 1) = a[0];     B(2, 2) = 0.0;

            Erg = cos(gamma) * unity() + (1.0 - cos(gamma)) * A + sin(gamma) * B;
            return Erg;
        }

        Matrix <double> rotMatrixA(Vector<double> n, Vector<double> k, double gamma)
        {
            Vector<double> e0, e1, e2;
            Matrix<double> H, R, D, Erg;
            e2 = n;
            e0 = k;
            e1 = e2 % e0;
            trafo(e0, e1, e2, H, R);
            D = Dz(gamma);
            Erg = R * D * H;
            return Erg;
        }

        Matrix<double> Dz(double phi)
        {
            Matrix<double> Erg;
            Erg(0, 0) = cos(phi);  Erg(0, 1) = -sin(phi); Erg(0, 2) = 0.0;
            Erg(1, 0) = sin(phi);  Erg(1, 1) = cos(phi); Erg(1, 2) = 0.0;
            Erg(2, 0) = 0.0;       Erg(2, 1) = 0.0;      Erg(2, 2) = 1.0;
            return Erg;
        }

        Matrix<double> Dx(double phi)
        {
            Matrix<double> Erg;
            Erg(0, 0) = 1.0;
            Erg(1, 1) = cos(phi);  Erg(1, 2) = sin(phi);
            Erg(2, 1) = -sin(phi); Erg(2, 2) = cos(phi);
            return Erg;
        }

        Matrix<double> Dy(double phi)
        {
            /*
              Dreht um die y-Achse von der positiven z-Achse zur positiven x-Achse
            */

            Matrix<double> Erg;
            Erg(0, 0) = cos(phi);  Erg(0, 1) = 0.0; Erg(0, 2) = sin(phi);
            Erg(1, 0) = 0.0;      Erg(1, 1) = 1.0; Erg(1, 2) = 0.0;
            Erg(2, 0) = -sin(phi); Erg(2, 1) = 0.0; Erg(2, 2) = cos(phi);
            return Erg;
        }

        void trafo(const Vector<double>& e0,
            const Vector<double>& e1,
            const Vector<double>& e2,
            Matrix<double>& H,
            Matrix<double>& R)
            /*
              Berechnung von (H)in- und (R)ücktransformationsmatrix aus dem kartesischen
              Koordinatensystem in das von e0,e1 und e2 gebildete System

            */
        {
            /* H(0,0)=e0[0]; H(0,1)=e1[0]; H(0,2)=e2[0];
             H(1,0)=e0[1]; H(1,1)=e1[1]; H(1,2)=e2[1];
             H(2,0)=e0[2]; H(2,1)=e1[2]; H(2,2)=e2[2];
            */
            bool found;
            H(0, 0) = e0[0]; H(0, 1) = e0[1]; H(0, 2) = e0[2];
            H(1, 0) = e1[0]; H(1, 1) = e1[1]; H(1, 2) = e1[2];
            H(2, 0) = e2[0]; H(2, 1) = e2[1]; H(2, 2) = e2[2];

            R = invert(H, found);
        }

        Matrix<double> rotMatrix(Vector<double> P, double dtheta, double dphi)
        {
            Matrix<double> Erg;
            double phi;
            double cphi, cdphi, sphi, sdphi;
            double cdtheta, sdtheta;

            if (P[0] != 0.0) phi = atan2(P[1], P[0]); else phi = 0.0;
            if (phi < 0.0) phi += 2.0 * M_PI;

            cphi = cos(phi); cdphi = cos(dphi); sphi = sin(phi); sdphi = sin(dphi);
            cdtheta = cos(dtheta); sdtheta = sin(dtheta);

            Erg(0, 0) = cdphi * cdtheta;
            Erg(0, 1) = -sdphi * cdtheta;
            Erg(0, 2) = cphi * cdphi * sdtheta - sphi * sdphi * sdtheta;

            Erg(1, 0) = sdphi * cdtheta;
            Erg(1, 1) = cdphi * cdtheta;
            Erg(1, 2) = sphi * cdphi * sdtheta + cphi * sdphi * sdtheta;

            Erg(2, 0) = -sdtheta / cphi;
            Erg(2, 1) = 0.0;
            Erg(2, 2) = cdtheta;

            return Erg;
        }

        Matrix<double> rotMatrix(double alpha, double beta, double gamma)
        {
            Matrix<double> M;
            double ca, cb, cg;
            double sa, sb, sg;
            ca = cos(alpha); cb = cos(beta); cg = cos(gamma);
            sa = sin(alpha); sb = sin(beta); sg = sin(gamma);

            M(0, 0) = cg * cb;   M(0, 1) = sa * sb * cg + ca * sg; M(0, 2) = -ca * sb * cg + sa * sg;
            M(1, 0) = -cb * sg;   M(1, 1) = -sa * sb * sg + ca * cg; M(1, 2) = ca * sb * sg + cg * sa;
            M(2, 0) = sb;   M(2, 1) = -sa * cb;          M(2, 2) = ca * cb;
            return M;
        }
    }
}