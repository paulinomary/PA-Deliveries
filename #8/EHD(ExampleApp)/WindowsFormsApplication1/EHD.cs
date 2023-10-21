using System;
using System.Collections.Generic;
using System.Text;
using System.Drawing;
using System.Drawing.Imaging;

namespace WindowsFormsApplication1
{
    class EHD
    {

       
        private int Thresshold;

        public EHD( int Thresshold)
        {
            this.Thresshold = Thresshold;
        }

        public double[] Apply(Bitmap srcImg)
        {
            double[] EHDTable = new double[80];
            int width = srcImg.Width;
            int height = srcImg.Height;
            double[,] ImageGrid = new double[width, height];

            for (int R = 0; R < 80; R++)
            {
                EHDTable[R] = 0;

            }

            PixelFormat fmt = (srcImg.PixelFormat == PixelFormat.Format8bppIndexed) ?
                       PixelFormat.Format8bppIndexed : PixelFormat.Format24bppRgb;

            BitmapData srcData = srcImg.LockBits(
               new Rectangle(0, 0, width, height),
               ImageLockMode.ReadOnly, fmt);



            int offset = srcData.Stride - ((fmt == PixelFormat.Format8bppIndexed) ? width : width * 3);

            unsafe
            {
                byte* src = (byte*)srcData.Scan0.ToPointer();


                for (int y = 0; y < height; y++)
                {
                    for (int x = 0; x < width; x++, src += 3)
                    {


                        int mean = (int)(0.114 * src[0] + 0.587 * src[1] + 0.299 * src[2]);


                        ImageGrid[x, y] = mean;


                    }

                    src += offset;


                }

            }

            srcImg.UnlockBits(srcData);

            //παραγωγη EDH
            double F1, F2, F3, F4, F5 = 0;
            int Possition = 0;
            double Maximum = 0;

            int W =Convert.ToInt32((double) width / 4);
            int H = Convert.ToInt32((double)height / 4);

            for (int i = 0; i < width - 2; i += 2)
            {
                for (int j = 0; j < height - 2; j += 2)
                {
                    F1 = F2 = F3 = F4 = F5 = 0;

                    F1 = ImageGrid[i, j] * 1 + ImageGrid[i + 1, j] * -1 + ImageGrid[i, j + 1] * 1 + ImageGrid[i + 1, j + 1] * -1;
                    F2 = ImageGrid[i, j] * 1 + ImageGrid[i + 1, j] * 1 + ImageGrid[i, j + 1] * -1 + ImageGrid[i + 1, j + 1] * -1;
                    F3 = ImageGrid[i, j] * Math.Sqrt(2) + ImageGrid[i + 1, j] * 0 + ImageGrid[i, j + 1] * 0 + ImageGrid[i + 1, j + 1] * -Math.Sqrt(2);
                    F4 = ImageGrid[i, j] * 0 + ImageGrid[i + 1, j] * Math.Sqrt(2) + ImageGrid[i, j + 1] * -Math.Sqrt(2) + ImageGrid[i + 1, j + 1] * 0;
                    F5 = ImageGrid[i, j] * 2 + ImageGrid[i + 1, j] * -2 + ImageGrid[i, j + 1] * -2 + ImageGrid[i + 1, j + 1] * 2;

                    Maximum = Math.Max(F1, Math.Max(F2, Math.Max(F3, Math.Max(F4, F5))));

                    if (Maximum == F1) Possition = 0;
                    if (Maximum == F2) Possition = 1;
                    if (Maximum == F3) Possition = 2;
                    if (Maximum == F4) Possition = 3;
                    if (Maximum == F5) Possition = 4; 


                    if (Maximum >= Thresshold)
                    {

                        EHDTable[5 * ((int)(4 * (i / W)) + (int)(j / H)) + Possition]++;


                    }


                }

            }

            double NormSum = 0;

            for (int i = 0; i < 80; i++)
            {

                NormSum += EHDTable[i];
            }

            if (NormSum > 0)
            {
                for (int i = 0; i < 80; i++)
                {

                    EHDTable[i] = EHDTable[i] / NormSum;
                }
            }

            
            return (EHDTable);

        }
       



    }
}
