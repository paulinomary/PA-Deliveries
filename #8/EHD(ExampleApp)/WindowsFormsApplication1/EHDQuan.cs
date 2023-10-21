using System;
using System.Collections.Generic;
using System.Text;
using System.Drawing;
using System.Drawing.Imaging;

namespace WindowsFormsApplication1
{
    class EHDQuant
    {

        private static double[,] QuantTable =
                    {{0.010867, 0.057915, 0.099526, 0.144849, 0.195573, 0.260504, 0.358031, 0.530128},
                    {0.012266, 0.069934, 0.125879, 0.182307, 0.243396, 0.314563, 0.411728, 0.564319},
                    {0.004193, 0.025852, 0.046860, 0.068519, 0.093286, 0.123490, 0.161505, 0.228960},
                    {0.004174, 0.025924, 0.046232, 0.067163, 0.089655, 0.115391, 0.151904, 0.217745},
                    {0.006778, 0.051667, 0.108650, 0.166257, 0.224226, 0.285691, 0.356375, 0.450972}};



        public double[] Apply(double[] Local_Edge_Histogram)
        {
            double[] Edge_HistogramElement = new double[Local_Edge_Histogram.Length];
            double iQuantValue = 0;

            for (int i = 0; i < Local_Edge_Histogram.Length; i++)
            {
                for (int j = 0; j < 8; j++)
                {
                    Edge_HistogramElement[i] = j;
                    if (j < 7)
                        iQuantValue = (QuantTable[i % 5, j] + QuantTable[i % 5, j + 1]) / 2.0;
                    else
                        iQuantValue = 1.0;
                    if (Local_Edge_Histogram[i] <= iQuantValue)
                    {
                        break;
                    }
                }

            }
            return Edge_HistogramElement;
        }
    }
}
