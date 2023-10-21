
using System;
using System.Drawing;
using System.Drawing.Imaging;



/*
* This file is a modified part of the Caliph and Emir project: http://www.SemanticMetadata.net.
*
* Caliph & Emir is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 2 of the License, or
* (at your option) any later version.
*
* Caliph & Emir is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with Caliph & Emir; if not, write to the Free Software
* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/


/*
* DominantColorLogic.java
*
* Created on 25. June 2004, 07:15 and updated until 04. November 2004
* Author: Wolfgang Seiringer, 002521, w.seiringer@utanet.at for the vizir-project (visual information retrieval project)
* (http://vizir.ims.tuwien.ac.at) project of the IMS (Interactive Media Systems Group) www.ims.tuwien.ac.at, at the TU Wien
* www.tuwien.ac.at
 * 
 * 
 * 
 * * C# Version
 * This Code is a modification of Caliph& Emir Project.
 * Part of img(Rummager) project
 * © 2006-2008 Savvas Chatzichristofis
 * http://savvash.blogspot.com
 * savvash@gmail.com, schatzic@ee.duth.gr
 * If you use this code please cite:
 * Mathias Lux, S. A. Chatzichristofis, "LIRe: Lucene Image Retrieval - An Extensible Java CBIR Library", ACM International Conference on Multimedia 2008, Vancouver, BC, Canada October 27 – 31, 2008, Open Source Application Competition.

 * 
*/


public class DCD_Descriptor
{

    private int imageHeight;
    private int imageWidth;
    private int imsize;
    //--------------------------------------------------------//

    // Holding the actual image to process:
    private Bitmap Bitmap = null;

    private int SC_BIT = 5;
    //    private final char quantImageAlpha = 255;
    private static int FLT_MAX = 1000000;
    private float DSTMIN = (float)255.0;
    private float SPLFCT = (float)0.10;
    private int DCNMAX = 8;
    private double EPSGLOB = 0.01; //Global stopping criteria, see extract method
    private double EPSSPLT = 0.02; //Criteria for splitting the colorbins
    private double VARTHR = 50.0;
    private int m_Coherency; //saves the spatical coherency

    //    private double dist, distold = FLT_MAX, distnew, eps = 1.0;
    private int j;

    int m_MaxSize = 8;
    int m_CurrSize = 1;

    private float[] m_Weights = new float[8]; //Saves the Percentages of the Dominant Colors
    private float[][] m_Centroids = new float[8][];


    //Saves the values of the Domiant Colors
    //      private int[] closest = new int[imsize]; //Is used for Clustering, see Method Cluster
    private int[] closest; //Is used for Clustering, see Method Cluster

    private int component0, component1, component2;
    private int bin_number0, bin_number1, bin_number2;

    // Variables for Saving the Values extracted from the XML-String
    public int DominantSize = 0; //Saves the "Size" extracted from the XML-String
    public int spatcoher = 0; //Saves the "SpatialCoherency" extracted from the XML-String
    public int[] Percentage = null; //Saves the "Percentage" extracted from the XML-String
    public int[][] ColorValueIndex = null; //Saves the "ColorValueIndex" extracted from the XML-String
    //-----------------

    //----- Variables used for the dominant color parameters
    private int VariancePresent = 0;
    private int SpatialCoherency = 0;
    private int ColorSpacePresentInt = 0;
    private String ColorSpacePresentStr = "0";
    private int ColorQuantizationPresent = 0;


    private int[][] icnts_f;
    private int[] iwgts_f;
    public float[][] m_Variances = new float[8][];
    public int[] interleavedArray;

    //---------------------------------------------------//

    public int[] percentages;
    public int[][] colorValues;


    public void extractDescriptor(Bitmap image)
    {
        this.imageHeight = image.Height;
        this.imageWidth = image.Width;
        this.imsize = this.imageHeight * this.imageWidth;
        this.closest = new int[imsize];
        this.Bitmap = image;

        this.extractFeature();

        for (int i = 0; i < m_CurrSize; i++)
        {
            colorValues[i] = new int[3];

            percentages[i] = iwgts_f[i];
            colorValues[i][0] = icnts_f[i][0];
            colorValues[i][1] = icnts_f[i][1];
            colorValues[i][2] = icnts_f[i][2];
        }


    }

    public void extractFeature()
    {

        int[] RGB;
        float[] LUV;

        if (Bitmap != null)
        {
            RGB = interleave(imsize, this.imageWidth, this.imageHeight);
        }
        else
        {
            RGB = interleavedArray;
        }

        LUV = new float[3 * imsize];

        rgb2luv(RGB, LUV, 3 * imsize); //Calculate the needed LUV-Values for the further operations

        char quantImageAlpha = Convert.ToChar(255);
        Extract(LUV, imsize, quantImageAlpha);//Calls the Main Extracting Method

    }//End Method extractFeature


    private static void rgb2luv(int[] RGB, float[] LUV, int size)
    {

        int i;
        double x, y, X, Y, Z, den, u2, v2, X0, Z0, Y0, u20, v20, r, g, b;

        X0 = (0.607 + 0.174 + 0.201);
        Y0 = (0.299 + 0.587 + 0.114);
        Z0 = (0.066 + 1.117);

        u20 = 4d * X0 / (X0 + 15d * Y0 + 3d * Z0);
        v20 = 9d * Y0 / (X0 + 15d * Y0 + 3d * Z0);

        for (i = 0; i < size; i += 3)
        {
            if (RGB[i] <= 20)
                r = (double)(8.715e-4 * RGB[i]);
            else
                r = (double)Math.Pow((RGB[i] + 25.245) / 280.245, 2.22);

            if (RGB[i + 1] <= 20)
                g = (double)(8.715e-4 * RGB[i + 1]);
            else
                g = (double)Math.Pow((RGB[i + 1] + 25.245) / 280.245, 2.22);

            if (RGB[i + 2] <= 20)
                b = (double)(8.715e-4 * RGB[i + 2]);
            else
                b = (double)Math.Pow((RGB[i + 2] + 25.245) / 280.245, 2.22);

            X = 0.412453 * r + 0.357580 * g + 0.180423 * b;
            Y = 0.212671 * r + 0.715160 * g + 0.072169 * b;
            Z = 0.019334 * r + 0.119193 * g + 0.950227 * b;

            if (X == 0.0 && Y == 0.0 && Z == 0.0)
            {
                x = 1.0 / 3.0;
                y = 1.0 / 3.0;
            }
            else
            {
                den = X + Y + Z;
                x = X / den;
                y = Y / den;
            }

            den = -2 * x + 12 * y + 3;
            u2 = 4 * x / den;
            v2 = 9 * y / den;

            if (Y > 0.008856)
                LUV[i] = (float)(116 * Math.Pow(Y, 1.0 / 3.0) - 16);
            else
                LUV[i] = (float)(903.3 * Y);
            LUV[i + 1] = (float)(13 * LUV[i] * (u2 - u20));
            LUV[i + 2] = (float)(13 * LUV[i] * (v2 - v20));
        }
    }// End Method RGBtoLUV

    // --------------------- LUV to RGB conversion------------------

    /**
     * Converts the LUV-Values in the corresponding RGB Values (Range 256)
     */
    private static void luv2rgb(int[] RGB, float[] LUV, int size)
    {
        int i, k;
        double x, y, X, Y, Z, den, u2, v2, X0, Z0, Y0, u20, v20;
        double[] vec = new double[3];

        X0 = (0.607 + 0.174 + 0.201);
        Y0 = (0.299 + 0.587 + 0.114);
        Z0 = (0.066 + 1.117);


        u20 = 4 * X0 / (X0 + 15 * Y0 + 3 * Z0);
        v20 = 9 * Y0 / (X0 + 15 * Y0 + 3 * Z0);

        for (i = 0; i < size; i += 3)
        {
            if (LUV[i] > 0)
            {
                if (LUV[i] < 8.0)
                    Y = ((double)LUV[i]) / 903.3;
                else
                    Y = Math.Pow((((double)LUV[i]) + 16) / 116.0, 3.0);
                u2 = ((double)LUV[i + 1]) / 13.0 / ((double)LUV[i]) + u20;
                v2 = ((double)LUV[i + 2]) / 13.0 / ((double)LUV[i]) + v20;

                den = 6 + 3 * u2 - 8 * v2;
                x = 4.5 * u2 / den;
                y = 2.0 * v2 / den;

                X = x / y * Y;
                Z = (1 - x - y) / y * Y;
            }
            else
            {
                X = 0.0;
                Y = 0.0;
                Z = 0.0;
            }

            vec[0] = (3.240479 * X - 1.537150 * Y - 0.498536 * Z);
            vec[1] = (-0.969256 * X + 1.875992 * Y + 0.041556 * Z);
            vec[2] = (0.055648 * X - 0.204043 * Y + 1.057311 * Z);
            for (k = 0; k < 3; k++)
            {
                if (vec[k] <= 0.018)
                    vec[k] = 255 * 4.5 * vec[k];
                else
                    vec[k] = 255 * (1.099 * Math.Pow(vec[k], 0.45) - 0.099);
                if (vec[k] > 255)
                    vec[k] = 255;
                else if (vec[k] < 0) vec[k] = 0;
                RGB[i + k] = (int)Math.Round((float)vec[k]);
            }
        }
    }


    private void Extract(float[] imdata, int imsize, char quantImageAlpha)
    {

        double dist, distold = FLT_MAX, distnew, eps = 1.0, tmp;
        float aglfct = DSTMIN;
        float splfct = SPLFCT;
        int i, k;

        m_Centroids[0] = new float[3];
        m_Centroids[1] = new float[3];
        m_Centroids[2] = new float[3];
        m_Centroids[3] = new float[3];
        m_Centroids[4] = new float[3];
        m_Centroids[5] = new float[3];
        m_Centroids[6] = new float[3];
        m_Centroids[7] = new float[3];

        m_MaxSize = DCNMAX;
        m_CurrSize = 1;

        for (i = 0; i < m_MaxSize; i++)
        {
            m_Variances[i] = new float[3];


            m_Weights[i] = (float)0.0;

            m_Centroids[i][0] = (float)0.0;
            m_Centroids[i][1] = (float)0.0;
            m_Centroids[i][2] = (float)0.0;

            m_Variances[i][0] = (float)0.0;
            m_Variances[i][1] = (float)0.0;
            m_Variances[i][2] = (float)0.0;
        }


        i = 0;

        distnew = Cluster(closest, imdata, imsize, quantImageAlpha);

        while (eps > EPSGLOB)
        {

            //find centroids
            Centroids(closest, imdata, imsize, quantImageAlpha);

            //classify bins
            distnew = Cluster(closest, imdata, imsize, quantImageAlpha);

            // calculate total distortion
            if (distold > 0.0)
                eps = (distold - distnew) / distold;
            else
                eps = 0.0;
            distold = distnew;

            //decide on splitting
            if (i == 0 || ((eps < EPSSPLT) && (m_CurrSize < m_MaxSize)))
            {

                Split(closest, imdata, imsize, splfct, quantImageAlpha);
                distnew = Cluster(closest, imdata, imsize, quantImageAlpha);
                eps = 1.0;
            }

            // check for identical codevectors
            for (j = 0; j < m_CurrSize; j++)
            {
                for (k = 0; k < j; k++)
                {
                    dist = 0.0;
                    tmp = m_Centroids[j][0] - m_Centroids[k][0];
                    dist += tmp * tmp;
                    tmp = m_Centroids[j][1] - m_Centroids[k][1];
                    dist += tmp * tmp;
                    tmp = m_Centroids[j][2] - m_Centroids[k][2];
                    dist += tmp * tmp;
                }
            }

            i++;

        }
        //  System.out.println("Extract: iterations finished " + i);

        // merging using agglomerative clustering
        Agglom(aglfct);

        // calculate variances and normalise
        distnew = Cluster(closest, imdata, imsize, quantImageAlpha);
        Centroids(closest, imdata, imsize, quantImageAlpha);
        distnew = Cluster(closest, imdata, imsize, quantImageAlpha);

        if (this.VariancePresent != 0)
        {
            Vars(closest, imdata, imsize, quantImageAlpha);
        }

        for (j = 0; j < m_CurrSize; j++)
            m_Weights[j] /= imsize;

        // quantise and set descriptor members
        int[] iwgts = new int[m_CurrSize];


        int[][] icnts = new int[3 * m_CurrSize][]; ;
        for (i = 0; i < m_CurrSize; i++)
            icnts[i] = new int[3];

        int[][] ivars = new int[3 * m_CurrSize][];

        for (i = 0; i < m_CurrSize; i++)
            ivars[i] = new int[3];

        for (i = 0; i < m_CurrSize; i++)
        {
            iwgts[i] = (int)(31.9999 * m_Weights[i]);
            luv2rgb(icnts[i], m_Centroids[i], 3);
        }

        for (i = 0; i < m_CurrSize; i++)
        {
            for (int j = 0; j < 3; j++)
                ivars[i][j] = (m_Variances[i][j] > VARTHR) ? 1 : 0;
        }

        // calculate spatial coherency
        if (this.SpatialCoherency != 0)
        {
            m_Coherency = GetSpatialCoherency(imdata, 3, m_CurrSize, m_Centroids, quantImageAlpha);
        }
        else
            m_Coherency = 0;

        String colspcpres = this.ColorSpacePresentStr;
        //        int colspcpresInt = this.ColorSpacePresentInt;
        int colquantpres = this.ColorQuantizationPresent;

        int max_h;
        //        float conversionmatrix;

        component0 = 0;
        component1 = 1;
        component2 = 2; // default RGB

        bin_number0 = 32;
        bin_number1 = 32;
        bin_number2 = 32; // default 5 bit

        max_h = 256;

        if (colspcpres != "0")
        {

            int[] inI = new int[25];

            //            int[] res = in;

            for (i = 0; i < m_CurrSize; i++)
            {
                inI[12] = icnts[i][1];
                inI[12] = icnts[i][2];
                inI[12] = icnts[i][0];
            }
        }

        if (colquantpres != 0)
        {
            bin_number0 = this.GetBinNumberByComponent(component0);
            bin_number1 = this.GetBinNumberByComponent(component1);
            bin_number2 = this.GetBinNumberByComponent(component2);
        }

        for (i = 0; i < m_CurrSize; i++)
        {
            icnts[i][0] = icnts[i][0] * bin_number0 / max_h;
            icnts[i][1] = icnts[i][1] * bin_number1 >> 8;
            icnts[i][2] = icnts[i][2] * bin_number2 >> 8;
        }


        for (int we = 0; we < m_CurrSize; we++)
        {
            //      System.out.println("Percentage: " + iwgts[we]);
            //      System.out.println("Values: " + icnts[we][0] + " " + icnts[we][1] + " " + icnts[we][2]);
            //      System.out.println("------------------------------------------------");
        }
        //      ********************************************************************

        icnts_f = icnts;
        iwgts_f = iwgts;

    } //end method extract

    //---------------------------------------------------------------------------

    /**
     * Each color is asign a cluster
     */
    private double Cluster(int[] closest, float[] imdata, int imsize, char quantImageAlpha)
    {

        int i, j, jmin;
        double d1, d2, d3, dist, distmin, disttot = 0.0;
        float[] im1 = imdata;
        //        float[] im2 = imdata, im3 = imdata;
        int imsize_msk;
        //        char pAlpha;

        // now we cluster
        imsize_msk = 0;

        for (i = 0; i < imsize; i++)
        {
            jmin = 0;
            distmin = FLT_MAX;

            for (j = 0; j < m_CurrSize; j++)
            {
                d1 = im1[3 * i] - m_Centroids[j][0];
                d2 = im1[3 * i + 1] - m_Centroids[j][1];
                d3 = im1[3 * i + 2] - m_Centroids[j][2];

                dist = d1 * d1 + d2 * d2 + d3 * d3;
                if (dist < distmin)
                {
                    jmin = j;
                    distmin = dist;
                }
            }

            closest[i] = jmin;
            if (closest[i] != 0) disttot += distmin;
            imsize_msk++;

        }
        return disttot / imsize_msk;
    } // end method cluster
    //----------------------------------Ende Cluster-----------------------------------------

    private void Centroids(int[] closest, float[] imdata, int imsize, char quantImageAlpha)
    {
        int i, j;
        double weight;
        float[] im1 = imdata;
        //        float[] im2 = imdata, im3 = imdata;
        //        char pAlpha;

        //set the centroids values to 0
        for (j = 0; j < m_CurrSize; j++)
        {
            m_Weights[j] = (float)0.0;
            m_Centroids[j][0] = (float)0.0;
            m_Centroids[j][1] = (float)0.0;
            m_Centroids[j][2] = (float)0.0;
        }

        for (i = 0; i < imsize; i++)
        {
            int ii = closest[i];
            m_Weights[ii]++;

            m_Centroids[ii][0] += im1[3 * i];
            m_Centroids[ii][1] += im1[3 * i + 1];
            m_Centroids[ii][2] += im1[3 * i + 2];
        }
        for (j = 0; j < m_CurrSize; j++)
        {
            weight = m_Weights[j];

            if (weight != 0.0)
            {
                m_Centroids[j][0] /= (float)weight;
                m_Centroids[j][1] /= (float)weight;
                m_Centroids[j][2] /= (float)weight;
            }
        }

    } //end method centroids
    //------------------------- END CENTROIDS --------------------------------------------------

    private void Vars(int[] closest, float[] imdata, int imsize, char quantImageAlpha)
    {
        int i, j;
        double tmp;

        //        char pAlpha;

        //set the m_Variances Values to 0
        for (i = 0; i < m_CurrSize; i++)
        {
            m_Variances[i][0] = (float)0.0;
            m_Variances[i][1] = (float)0.0;
            m_Variances[i][2] = (float)0.0;
        }

        for (i = 0; i < imsize; i++)
        {
            j = closest[i];
            tmp = imdata[3 * i] - m_Centroids[j][0];
            m_Variances[j][0] += (float)(tmp * tmp);
            tmp = imdata[3 * i + 1] - m_Centroids[j][1];
            m_Variances[j][1] += (float)(tmp * tmp);
            tmp = imdata[3 * i + 2] - m_Centroids[j][2];
            m_Variances[j][2] += (float)(tmp * tmp);
        }

        //Normalise values
        for (j = 0; j < m_CurrSize; j++)
        {
            m_Variances[j][0] /= m_Weights[j];
            m_Variances[j][1] /= m_Weights[j];
            m_Variances[j][2] /= m_Weights[j];
        }

    } // end method vars

    private void Split(int[] closest, float[] imdata, int imsize, double factor, char quantImageAlpha)
    {

        int i, j, jmax = 0;
        double d1, d2, d3, diff1, diff2, diff3;

        double[] d1s = new double[8];
        double[] d2s = new double[8];
        double[] d3s = new double[8];
        double[] dists = new double[8];
        double distmax = 0.0;

        //        char pAlpha;

        //set distortion values to 0
        for (j = 0; j < m_CurrSize; j++)
        {
            d1s[j] = 0.0;
            d2s[j] = 0.0;
            d3s[j] = 0.0;
            dists[j] = 0.0;
        }

        for (i = 0; i < imsize; i++)
        {
            j = closest[i];
            d1 = imdata[3 * i] - m_Centroids[j][0];
            d2 = imdata[3 * i + 1] - m_Centroids[j][1];
            d3 = imdata[3 * i + 2] - m_Centroids[j][2];
            d1s[j] += d1 * d1;
            d2s[j] += d2 * d2;
            d3s[j] += d3 * d3;
        }

        for (j = 0; j < m_CurrSize; j++)
        {
            dists[j] = d1s[j] + d2s[j] + d3s[j];
            d1s[j] /= m_Weights[j];
            d2s[j] /= m_Weights[j];
            d3s[j] /= m_Weights[j];
        }

        for (j = 0; j < m_CurrSize; j++)
            if (dists[j] > distmax)
            {
                jmax = j;
                distmax = dists[j];
            }

        diff1 = factor * Math.Sqrt(d1s[jmax]);
        diff2 = factor * Math.Sqrt(d2s[jmax]);
        diff3 = factor * Math.Sqrt(d3s[jmax]);

        m_Centroids[m_CurrSize][0] = (float)(m_Centroids[jmax][0] + diff1);
        m_Centroids[m_CurrSize][1] = (float)(m_Centroids[jmax][1] + diff2);
        m_Centroids[m_CurrSize][2] = (float)(m_Centroids[jmax][2] + diff3);

        m_Centroids[jmax][0] = (float)(m_Centroids[jmax][0] - diff1);
        m_Centroids[jmax][1] = (float)(m_Centroids[jmax][1] - diff2);
        m_Centroids[jmax][2] = (float)(m_Centroids[jmax][2] - diff3);

        m_CurrSize++;

        percentages = new int[m_CurrSize];
        colorValues = new int[m_CurrSize][];
    } // end method Split


    private void Agglom(double distthr)
    {

        double d1, d2, d3, distmin = 0.0;
        double[][] dists = new double[8][];
        double w1min, w2min;
        int ja, jb = 0, jamin, jbmin;

        for (int i = 0; i < 8; i++)
        {
            dists[i] = new double[8];
        }


        do
        {

            for (ja = 0; ja < m_CurrSize; ja++)
                for (jb = 0; jb < ja; jb++)
                {
                    d1 = m_Centroids[ja][0] - m_Centroids[jb][0];
                    d2 = m_Centroids[ja][1] - m_Centroids[jb][1];
                    d3 = m_Centroids[ja][2] - m_Centroids[jb][2];
                    //                    dists[ja][jb] = d1 * d1 + d2 * d2 + d3 * d3;

                    /* START: added by janine to circumvent divisions by 0 */
                    //Distanz > 0, nur wenn weight nicht 0 ist!
                    if (m_Weights[ja] > 0.0 && m_Weights[jb] > 0.0)
                        dists[ja][jb] = (d1 * d1 + d2 * d2 + d3 * d3);
                    else
                        dists[ja][jb] = 0.0;
                    /* END: added by janine to circumvent divisions by 0 */

                }

            distmin = FLT_MAX;
            jamin = 0;
            jbmin = 0;
            for (ja = 0; ja < m_CurrSize; ja++)
                for (jb = 0; jb < ja; jb++)
                    if (dists[ja][jb] < distmin)
                    {
                        distmin = dists[ja][jb];
                        jamin = ja;
                        jbmin = jb;
                    }

            if (distmin > distthr)
                break;


            w1min = m_Weights[jamin];
            w2min = m_Weights[jbmin];

            // 'if' eingefuegt von Janine (division durch 0)
            if (w1min + w2min > 0)
            {
                m_Centroids[jbmin][0] = (float)((w1min * m_Centroids[jamin][0] + w2min * m_Centroids[jbmin][0]) / (w1min + w2min));
                m_Centroids[jbmin][1] = (float)((w1min * m_Centroids[jamin][1] + w2min * m_Centroids[jbmin][1]) / (w1min + w2min));
                m_Centroids[jbmin][2] = (float)((w1min * m_Centroids[jamin][2] + w2min * m_Centroids[jbmin][2]) / (w1min + w2min));
            }
            m_Weights[jbmin] += (float)w1min;
            m_CurrSize--;


            for (jb = jamin; jb < m_CurrSize; jb++)
            {
                m_Weights[jb] = m_Weights[jb + 1];
                m_Centroids[jb][0] = m_Centroids[jb + 1][0];
                m_Centroids[jb][1] = m_Centroids[jb + 1][1];
                m_Centroids[jb][2] = m_Centroids[jb + 1][2];
                m_Variances[jb][0] = m_Variances[jb + 1][0];
                m_Variances[jb][1] = m_Variances[jb + 1][1];
                m_Variances[jb][2] = m_Variances[jb + 1][2];
            }

        } while (m_CurrSize > 1 && distmin < distthr);

    } // end method aglom


    private int GetSpatialCoherency(float[] ColorData, int dim, int N, float[][] col_float, char quantImageAlpha)
    {

        //        char pAlpha;
        double CM = 0.0;
        int NeighborRange = 1;

        float SimColorAllow = (float)Math.Sqrt(DSTMIN);
        Boolean[] IVisit = new Boolean[imsize];

        for (int x = 0; x < (imsize); x++)
        {
            IVisit[x] = false;
        }

        int All_Pixels = 0;

        {

            for (int x = 0; x < (imsize); x++)
                if (IVisit[x] == false) All_Pixels++;

        }

        {
            for (int i = 0; i < N; i++)
            {

                int Corres_Pixels = 0;
                double Coherency = GetCoherencyWithColorAllow(ColorData, dim, IVisit, col_float[i][0], col_float[i][1], col_float[i][2],
                        SimColorAllow, NeighborRange, Corres_Pixels++);

                CM += (Coherency * (double)Corres_Pixels) / 2.7;

            }
        }
        return QuantizeSC(CM);
    }//end method get spatial coherency

    private double GetCoherencyWithColorAllow(float[] ColorData, int dim, Boolean[] IVisit,

                                              float l, float u, float v,
                                              float Allow, int NeighborRange, int OUTPUT_Corres_Pixel_Count)
    {
        int count, i, j;
        int Neighbor_Count = 0;

        int Pixel_Count = 0;
        double Coherency = 0.0;


        int width = imageWidth;


        int height = imageHeight;
        int ISize = imsize * dim;

        for (count = 0; count < ISize; count += dim)
        {
            i = (count / dim) % width; //width
            j = (count / dim) / width; //height

            float l1, u1, v1;
            l1 = ColorData[count];
            u1 = ColorData[count + 1];
            v1 = ColorData[count + 2];


            double distance;

            distance = Math.Sqrt(Math.Pow(l - l1, 2) + Math.Pow(u - u1, 2) + Math.Pow(v - v1, 2));

            if ((distance < Allow) && (IVisit[count / dim] == false))//no overlap checking
            {
                IVisit[count / dim] = true;
                Pixel_Count++;
                int nSameNeighbor = 0;
                for (int y = j - NeighborRange; y <= j + NeighborRange; y++)
                {
                    for (int x = i - NeighborRange; x <= i + NeighborRange; x++)
                    {
                        if (!((i == x) && (j == y)))
                        {
                            int Index = (y * width + x) * dim;
                            if ((Index >= 0) && (Index < ISize))
                            {
                                float l2 = ColorData[Index];
                                float u2 = ColorData[Index + 1];
                                float v2 = ColorData[Index + 2];
                                distance = Math.Sqrt(Math.Pow((l - l2), 2) + Math.Pow((u - u2), 2) + Math.Pow((v - v2), 2));
                                if (distance < Allow)
                                {
                                    nSameNeighbor++;
                                }
                            }
                        }
                    }
                }
                Neighbor_Count += nSameNeighbor;
            }
        }

        OUTPUT_Corres_Pixel_Count = Pixel_Count;

        int neighbor_check_window_size = (NeighborRange << 1) + 1;
        neighbor_check_window_size *= neighbor_check_window_size;

        if (Pixel_Count == 0)
            Coherency = 0.0;
        else
            Coherency = (double)Neighbor_Count / (double)Pixel_Count / (double)(neighbor_check_window_size - 1);

        return Coherency;
    }// method get Coherency with Color Allow


    private int[] interleave(int size, int width, int height)
    {
        int[] pixelarray = new int[3 * size];
        int j = 0;

        ///


        PixelFormat fmt = (Bitmap.PixelFormat == PixelFormat.Format8bppIndexed) ?
                         PixelFormat.Format8bppIndexed : PixelFormat.Format24bppRgb;

        BitmapData srcData = Bitmap.LockBits(
           new Rectangle(0, 0, width, height),
           ImageLockMode.ReadOnly, fmt);

        int offset = srcData.Stride - srcData.Width * 3;


        unsafe
        {
            byte* src = (byte*)srcData.Scan0.ToPointer();

            for (int y = 0; y < height; y++)
            {
                for (int x = 0; x < width; x++, src += 3)
                {


                    pixelarray[3 * j] = src[2];
                    pixelarray[3 * j + 1] = src[1];
                    pixelarray[3 * j + 2] = src[0];
                    j++;
                }

                src += offset;


            }

        }

        Bitmap.UnlockBits(srcData);


        return pixelarray;
    }

    private int QuantizeSC(double sc)
    {
        if (sc < 0.70)
            return 1;
        else
            return (int)((sc - 0.70) / (1.0 - 0.70) * (Math.Pow(2.0, (double)SC_BIT) - 3.0) + .5) + 2;
    }

    private int GetBinNumberByComponent(int component)
    {
        int i;
        int NumComponents = 3;
        int[] m_Component = new int[3];
        int[] m_BinNumber = new int[3];

        int ColorSpace = this.ColorSpacePresentInt;
        String ColorSpaceStr = this.ColorSpacePresentStr;

        if (ColorSpaceStr == "Monochrome") NumComponents = 1;         // monochrome (0101)
        if (ColorSpace == 5) NumComponents = 1;                                 // monochrome (0101)

        for (i = 0; i < NumComponents; i++)
            if (m_Component[i] == component) return m_BinNumber[i];

        return -1;
    }






}
