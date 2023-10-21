// This code is based on the Caliph and Emir project: http://www.SemanticMetadata.net. */
//(c) Mathias Lux under GLU using Java
// Edit by Savvas Chatzichristofis - (c) http://www.chatzichristofis.info - Using C#
    
/*
 *   This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 */


using System;
using System.Collections.Generic;
using System.Text;

namespace Similarity_Matching
{
    class MPEG_7Simillarity
    {
            private static double[,] QuantTable =
            {{0.010867, 0.057915, 0.099526, 0.144849, 0.195573, 0.260504, 0.358031, 0.530128},
                    {0.012266, 0.069934, 0.125879, 0.182307, 0.243396, 0.314563, 0.411728, 0.564319},
                    {0.004193, 0.025852, 0.046860, 0.068519, 0.093286, 0.123490, 0.161505, 0.228960},
                    {0.004174, 0.025924, 0.046232, 0.067163, 0.089655, 0.115391, 0.151904, 0.217745},
                    {0.006778, 0.051667, 0.108650, 0.166257, 0.224226, 0.285691, 0.356375, 0.450972}};
        protected static int[,] weightMatrix = new int[3, 64];

        public double calculateSCDDistance(double[] HistogramA, double[] HistogramB)
        {

            double result = 0;

            for (int l = 0; l < 64; l++)
            {
                result += Math.Abs(HistogramA[l] - HistogramB[l]);

            }
            return result;
        }

        public double calculateEHDDistance(double[] edgeHistogramA, double[] edgeHistogramB)
        {

            double result = 0;

            for (int i = 0; i < edgeHistogramA.Length; i++)
            {
                result += Math.Abs((float)QuantTable[i % 5, (int)edgeHistogramA[i]] - (float)QuantTable[i % 5, (int)edgeHistogramB[i]]);
            }
/*
            for (int i = 0; i <= 4; i++)
            {
                result += 5f * Math.Abs((float)edgeHistogramA[i] - (float)edgeHistogramB[i]);
            }
            for (int i = 5; i < 80; i++)
            {
                result += Math.Abs((float)edgeHistogramA[i] - (float)edgeHistogramB[i]);
            }
*/
            return result;
        }

        private static void setWeightingValues()
        {
            weightMatrix[0, 0] = 2;
            weightMatrix[0, 1] = weightMatrix[0, 2] = 2;
            weightMatrix[1, 0] = 2;
            weightMatrix[1, 1] = weightMatrix[1, 2] = 1;
            weightMatrix[2, 0] = 4;
            weightMatrix[2, 1] = weightMatrix[2, 2] = 2;
            for (int i = 0; i < 3; i++)
            {
                for (int j = 3; j < 64; j++)
                    weightMatrix[i, j] = 1;
            }
        }

        public  double calculateCLDDistance(int[] YCoeff1, int[] CbCoeff1, int[] CrCoeff1, int[] YCoeff2, int[] CbCoeff2, int[] CrCoeff2)
        {
            int numYCoeff1, numYCoeff2, CCoeff1, CCoeff2, YCoeff, CCoeff;
            //Numbers of the Coefficients of two descriptor values.
            numYCoeff1 = YCoeff1.Length;
            numYCoeff2 = YCoeff2.Length;
            CCoeff1 = CbCoeff1.Length;
            CCoeff2 = CbCoeff2.Length;
            //take the Minimal Coeff-number
            YCoeff = Math.Min(numYCoeff1, numYCoeff2);
            CCoeff = Math.Min(CCoeff1, CCoeff2);
            setWeightingValues();
            int j;
            int[] sum = new int[3];
            int diff;
            sum[0] = 0;
            for (j = 0; j < YCoeff; j++)
            {
                diff = (YCoeff1[j] - YCoeff2[j]);
                sum[0] += (weightMatrix[0, j] * diff * diff);
            }
            sum[1] = 0;
            for (j = 0; j < CCoeff; j++)
            {
                diff = (CbCoeff1[j] - CbCoeff2[j]);
                sum[1] += (weightMatrix[1, j] * diff * diff);
            }
            sum[2] = 0;
            for (j = 0; j < CCoeff; j++)
            {
                diff = (CrCoeff1[j] - CrCoeff2[j]);
                sum[2] += (weightMatrix[2, j] * diff * diff);
            }
            //returns the distance between the two desciptor values
            return Math.Sqrt(sum[0] * 1.0) + Math.Sqrt(sum[1] * 1.0) + Math.Sqrt(sum[2] * 1.0);
        }



 

    }
}
