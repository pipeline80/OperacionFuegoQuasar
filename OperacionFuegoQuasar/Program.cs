using MathNet.Numerics.LinearAlgebra.Complex;
using System.Numerics;
using System;

namespace OperacionFuegoQuasar
{
    class Program
    {
        private const int V = 200; private const int V1 = 200;

        static void Main(string[] args)
        {
            //The phone's x and y coordinates are pre-set



            //Point p1 = new Point(-500, -200, 100);
            //Point p2 = new Point(100, -100, 115.5);
            //Point p3 = new Point(500, 100, 142.7);
            //double[] a = Trilateration.Compute(p1, p2, p3);
            //Console.WriteLine("LatLon: " + a[0] + ", " + a[1]);

            //point p1 = new point(); p1.x = -500; p1.y = -200; p1.r = 100;
            //point p2 = new point(); p2.x = 100; p2.y = -100; p2.r = 115.5;
            //point p3 = new point(); p3.x = 500; p3.y = 100; p3.r = 142.7;

            point p1 = new point(); p1.x = -19.6685; p1.y = -69.1942; p1.r = 84;
            point p2 = new point(); p2.x = -20.2705; p2.y = -70.1311; p2.r = 114;
            point p3 = new point(); p3.x = -20.5656; p3.y = -70.1807; p3.r = 120;

            point a = Trilateration.Computes(p1, p2, p3);
            Console.WriteLine("LatLon: " + a.x + ", " + a.y);
        }

       
    }
    public class Point
    {
        private double lt, ln, r;
        public Point(double lt, double ln, double r) { this.lt = lt; this.ln = ln; this.r = r; }
        public double glt() { return lt; }
        public double gln() { return ln; }
        public double gr() { return r; }
    }

    public class point
    {
        public double x, y, r;
    }

    public class Trilateration
    {
        public static double norm(point p) // get the norm of a vector
        {
            return Math.Pow(Math.Pow(p.x, 2) + Math.Pow(p.y, 2), 0.5);
        }

        public static point Computes(point point1, point point2, point point3)
        {
            double p2p1Distance = Math.Pow(Math.Pow(point2.x - point1.x, 2) + Math.Pow(point2.y - point1.y, 2), 0.5);
            point ex = new point(); ex.x = (point2.x - point1.x) / p2p1Distance; ex.y = (point2.y - point1.y) / p2p1Distance;
            point aux = new point(); aux.x = point3.x - point1.x; aux.y = point3.y - point1.y;
            double i = ex.x * aux.x + ex.y * aux.y;
            point aux2 = new point(); aux2.x = point3.x - point1.x - i * ex.x; aux2.y = point3.y - point1.y - i * ex.y;
            point ey = new point(); ey.x = aux2.x / norm(aux2); ey.y = aux2.y / norm(aux2);

            double j = ey.x * aux.x + ey.y * aux.y;

            double x = (Math.Pow(point1.r, 2) - Math.Pow(point2.r, 2) + Math.Pow(p2p1Distance, 2)) / (2 * p2p1Distance);
            double y = (Math.Pow(point1.r, 2) - Math.Pow(point3.r, 2) + Math.Pow(i, 2) + Math.Pow(j, 2)) / (2 * j) - i * x / j;
            //result coordinates
            double finalX = point1.x + x * ex.x + y * ey.x;
            double finalY = point1.y + x * ex.y + y * ey.y;

            point returnPoint = new point();
            returnPoint.x = finalX; returnPoint.y = finalY;
            return returnPoint;
        }
        public static double[] Compute(Point p1, Point p2, Point p3)
        {
            double[] a = new double[3];
            double[] b = new double[3];
            double c, d, f, g, h;
            double[] i = new double[2];
            double k;
            c = p2.glt() - p1.glt();
            d = p2.gln() - p1.gln();
            f = (180 / Math.PI) * Math.Acos(Math.Abs(c) / Math.Abs(Math.Sqrt(Math.Pow(c, 2) + Math.Pow(d, 2))));
            if ((c > 0 && d > 0)) { f = 360 - f; }
            else if ((c < 0 && d > 0)) { f = 180 + f; }
            else if ((c < 0 && d < 0)) { f = 180 - f; }
            a = C(c, d, B(A(D(p2.gr()))), f);
            b = C(p3.glt() - p1.glt(), p3.gln() - p1.gln(), B(A(D(p3.gr()))), f);
            g = (Math.Pow(B(A(D(p1.gr()))), 2) - Math.Pow(a[2], 2) + Math.Pow(a[0], 2)) / (2 * a[0]);
            h = (Math.Pow(B(A(D(p1.gr()))), 2) - Math.Pow(b[2], 2) - Math.Pow(g, 2) + Math.Pow(g - b[0], 2) + Math.Pow(b[1], 2)) / (2 * b[1]);
            i = C(g, h, 0, -f);
            i[0] = i[0] + p1.glt() - 0.086;
            i[1] = i[1] + p1.gln() - 0.004;
            k = E(i[0], i[1], p1.glt(), p1.gln());
            if (k > p1.gr() * 2) { i = null; }
            else
            {
                if (i[0] < -90 || i[0] > 90 || i[1] < -180 || i[1] > 180) { i = null; }
            }
            return i;
        }
        private static double A(double a) { return a * 7.2; }
        private static double B(double a) { return a / 900000; }
        private static double[] C(double a, double b, double c, double d) { return new double[] { a * Math.Cos((Math.PI / 180) * d) - b * Math.Sin((Math.PI / 180) * d), a * Math.Sin((Math.PI / 180) * d) + b * Math.Cos((Math.PI / 180) * d), c }; }
        private static double D(double a) { return 730.24198315 + 52.33325511 * a + 1.35152407 * Math.Pow(a, 2) + 0.01481265 * Math.Pow(a, 3) + 0.00005900 * Math.Pow(a, 4) + 0.00541703 * 180; }
        private static double E(double a, double b, double c, double d) { double e = Math.PI, f = e * a / 180, g = e * c / 180, h = b - d, i = e * h / 180, j = Math.Sin(f) * Math.Sin(g) + Math.Cos(f) * Math.Cos(g) * Math.Cos(i); if (j > 1) { j = 1; } j = Math.Acos(j); j = j * 180 / e; j = j * 60 * 1.1515; j = j * 1.609344; return j; }
    }
}
