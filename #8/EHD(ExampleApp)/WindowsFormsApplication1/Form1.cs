using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Windows.Forms;
using System.Windows.Forms.DataVisualization.Charting;
using ZedGraph;

//Programed by George Constantinou georkons1@ee.duth.gr
//bugs report georkons1@ee.duth.gr



namespace WindowsFormsApplication1
{
    public partial class Form1 : Form
    {

        Bitmap Eikona;

        public Form1()
        {
            InitializeComponent();
        }

        private void Form1_Load(object sender, EventArgs e)
        {
            
        }

        private void button1_Click(object sender, EventArgs e)
        {
           
            openFileDialog1.ShowDialog();
            Eikona = new Bitmap(openFileDialog1.FileName);
            pictureBox1.Image = Eikona;
            groupBox1.Enabled = true;
           
        }

        private void FormatingBHDiagrams(ZedGraphControl PaneL, PointPairList CenterPoints, int Length)
        {
            BarItem Centers;


            PaneL.IsShowPointValues = true;
          
            PaneL.GraphPane.Legend.IsVisible = false;

            
           
            PaneL.GraphPane.CurveList.Clear();

            Centers = PaneL.GraphPane.AddBar("Fuzzy Luminosity Histogram", CenterPoints, Color.DarkGray);
           
            PaneL.AxisChange();
            PaneL.Invalidate();

        }


        private void chart1_Click(object sender, EventArgs e)
        {
            




        }

        private void button2_Click(object sender, EventArgs e)
        {
            PointPairList PointsOfOutput = new PointPairList();

            EHDQuant QuantizationEHD = new EHDQuant();
             

            EHD DescriptorEHD = new EHD(Convert.ToInt32(trackBar1.Value));
            double[] Histogram = DescriptorEHD.Apply(Eikona);

            double[] Histogram150;

            

            if (checkBox1.Checked)
            {

                EHD150 DescriptorEHD150 = new EHD150();
                Histogram150 = DescriptorEHD150.Apply(Histogram);

                if (checkBox2.Checked) Histogram150 = QuantizationEHD.Apply(Histogram150);

                for (int i = 0; i < Histogram150.Length; i++)
                {
                    PointsOfOutput.Add(i, Histogram150[i]);
                }


                FormatingBHDiagrams(zedGraphControl1, PointsOfOutput, Histogram150.Length);


            }
            else
            {

                if (checkBox2.Checked) Histogram = QuantizationEHD.Apply(Histogram);

                for (int i = 0; i < Histogram.Length; i++)
                {
                    PointsOfOutput.Add(i, Histogram[i]);
                }


                FormatingBHDiagrams(zedGraphControl1, PointsOfOutput, Histogram.Length);
            }
        }

        private void trackBar1_Scroll(object sender, EventArgs e)
        {
            label1.Text =trackBar1.Value.ToString();
        }

        private void button3_Click(object sender, EventArgs e)
        {
            
           
        }
    }
}
