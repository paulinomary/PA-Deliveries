using System;
using System.Collections.Generic;
using System.Text;
using System.Drawing;
using System.Drawing.Imaging;

namespace WindowsFormsApplication1
{
    class EHD150
    {

        public double[] Apply(double[] EHDHistogram)
        {
            double[] EHD150Diagram = new double[150];
            double[] F = new double[5];
            double TempSum = 0;
            double TempSum2 = 0;
            double TempSum3 = 0;
            double TempSum4 = 0;
            double TempSum5 = 0;

            //τα πρώτα 80 bins ειναι τα ίδια με το EHD
            for (int i = 0; i < 80; i++)
            {
                EHD150Diagram[i] = EHDHistogram[i];
            }

            // υπολογίζω τα bins κάθε κλίσης
            for (int i = 0; i < 5; i++)
            { // Αποφέυγω το διπλό for για θέμα χρόνου
                F[i] = EHDHistogram[0 + i] + EHDHistogram[5 + i] + EHDHistogram[10 + i] + EHDHistogram[15 + i] + EHDHistogram[20 + i] +
                       EHDHistogram[25 + i] + EHDHistogram[30 + i] + EHDHistogram[35 + i] + EHDHistogram[40 + i] + EHDHistogram[45 + i] +
                       EHDHistogram[50 + i] + EHDHistogram[55 + i] + EHDHistogram[60 + i] + EHDHistogram[65 + i] + EHDHistogram[70 + i] +
                       EHDHistogram[75 + i];

                EHD150Diagram[80 + i] = F[i];
                TempSum += F[i];
            }

            // Κονονικοποιώ τα Bins

            for (int i = 0; i < 5; i++)
            {
                EHD150Diagram[80 + i] = EHD150Diagram[80 + i] / TempSum;
            }



            // Υπολογίζω τα Semi Global Bins


            for (int i = 0; i < 20; i++) //Αρχικά υπολογίζω τις γραμμές και τις στίλες. Η αρίθμηση ακολουθεί αυτή του Paper
            { // Αποφέυγω το διπλό for για θέμα χρόνου
                EHD150Diagram[85 + i] = EHDHistogram[i] + EHDHistogram[20 + i] + EHDHistogram[40 + i] + EHDHistogram[60 + i];
                EHD150Diagram[105 + i] = EHDHistogram[i] + EHDHistogram[5 + i] + EHDHistogram[10 + i] + EHDHistogram[15 + i];
                TempSum2 += EHD150Diagram[85 + i];
                TempSum3 += EHD150Diagram[105 + i];
            }

            // Κονονικοποιώ τα Bins

            for (int i = 0; i < 20; i++)
            {
                EHD150Diagram[85 + i] = EHD150Diagram[85 + i] / (TempSum2 / 4);
                EHD150Diagram[105 + i] = EHD150Diagram[105 + i] / (TempSum3 / 4);
            }

            TempSum = 0;
            TempSum2 = 0;
            TempSum3 = 0;
            TempSum4 = 0;
            TempSum5 = 0;

            // 9 - 13
            for (int i = 0; i < 5; i++)
            { // Αποφέυγω το διπλό for για θέμα χρόνου
                EHD150Diagram[125 + i] = EHDHistogram[i] + EHDHistogram[5 + i] + EHDHistogram[20 + i] + EHDHistogram[25 + i];
                EHD150Diagram[130 + i] = EHDHistogram[10 + i] + EHDHistogram[15 + i] + EHDHistogram[30 + i] + EHDHistogram[35 + i];
                EHD150Diagram[135 + i] = EHDHistogram[40 + i] + EHDHistogram[45 + i] + EHDHistogram[60 + i] + EHDHistogram[65 + i];
                EHD150Diagram[140 + i] = EHDHistogram[50 + i] + EHDHistogram[55 + i] + EHDHistogram[70 + i] + EHDHistogram[75 + i];
                EHD150Diagram[145 + i] = EHDHistogram[25 + i] + EHDHistogram[30 + i] + EHDHistogram[45 + i] + EHDHistogram[50 + i];

                TempSum += EHD150Diagram[125 + i];
                TempSum2 += EHD150Diagram[130 + i];
                TempSum3 += EHD150Diagram[135 + i];
                TempSum4 += EHD150Diagram[140 + i];
                TempSum5 += EHD150Diagram[145 + i];

            }

            // Κονονικοποιώ τα Bins
            for (int i = 0; i < 5; i++)
            {
                EHD150Diagram[125 + i] = EHD150Diagram[125 + i] / TempSum;
                EHD150Diagram[130 + i] = EHD150Diagram[130 + i] / TempSum2;
                EHD150Diagram[135 + i] = EHD150Diagram[135 + i] / TempSum3;
                EHD150Diagram[140 + i] = EHD150Diagram[140 + i] / TempSum4;
                EHD150Diagram[145 + i] = EHD150Diagram[145 + i] / TempSum5;
            }






            return (EHD150Diagram);
        }



    }
}
