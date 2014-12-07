using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Windows.Forms;
using System.IO;
using System.Text.RegularExpressions;
using System.Drawing.Drawing2D;
using System.Drawing.Imaging;

namespace geo_demo
{
    public partial class Form1 : Form
    {
        enum DemoSet { Nunavut, Bible, USA }

        public Form1()
        {
            InitializeComponent();
        }

        private void Form1_Load(object sender, EventArgs e)
        {
            DemoSet demoSet = DemoSet.Nunavut;

            GenerateMapMasks(demoSet);
            GenerateHeatmaps(demoSet);

        }

        #region Generation of tiled images for display purposes.

        private void GenerateTiledImage(string outputDir, int numAcross, int numDown)
        {
            // Load tile file names

            string[] allFiles = Directory.GetFiles("output/" + outputDir);
            List<string> answerFiles = new List<string>();
            foreach (string file in allFiles)
            {
                if (file.EndsWith("_answer.png")) answerFiles.Add(file);
            }
            Bitmap example = new Bitmap(answerFiles[0]);
            int tileWidth = example.Width;
            int tileHeight = example.Height;
            int gutterWidth = 4;
            int gutterHeight = 70;


            // Generate tiled images
            
            int tilesPerImage = numAcross * numDown;
            int imageCount = (int)Math.Ceiling(answerFiles.Count() / (double)tilesPerImage);
            int tileIndex = 0;
            Font font = new Font("Arial", 28);
            StringFormat sf = new StringFormat(StringFormatFlags.NoClip);
            sf.LineAlignment = StringAlignment.Center;
            sf.Alignment = StringAlignment.Center;

            for (int i = 0; i < imageCount; i++)
            {
                Bitmap image = new Bitmap((tileWidth + gutterWidth) * numAcross - gutterWidth,
                    (tileHeight + gutterHeight) * numDown);
                Graphics g = Graphics.FromImage(image);
                g.Clear(Color.White);

                DrawTiledImage(numAcross, numDown, answerFiles, tileWidth, tileHeight, gutterWidth, gutterHeight, 
                    ref tileIndex, font, sf, g);

                image.Save("tiles/" + outputDir + "_" + i.ToString() + ".png", ImageFormat.Png);
                
            }

            
        }

        private static void DrawTiledImage(int numAcross, int numDown, List<string> answerFiles, int tileWidth, int tileHeight, int gutterWidth, int gutterHeight, ref int tileIndex, Font font, StringFormat sf, Graphics g)
        {
            for (int y = 0; y < numDown; y++)
            {
                for (int x = 0; x < numAcross; x++)
                {
                    if (tileIndex >= answerFiles.Count()) return;

                    // Draw tile
                    string tileFileName = answerFiles[tileIndex];
                    Bitmap tile = new Bitmap(tileFileName);
                    g.DrawImage(tile, x * (tileWidth + gutterWidth), y * (tileHeight + gutterHeight), tile.Width, tile.Height);

                    // Draw string
                    string toponym = Path.GetFileNameWithoutExtension(tileFileName).Replace("_answer", "");
                    g.DrawString(toponym, font, Brushes.Black, 
                        new RectangleF(x * (tileWidth + gutterWidth), (y * (tileHeight + gutterHeight)) + tileHeight + 2, tileWidth, gutterHeight/2), sf);

                    tileIndex++;
                }
            }
        }

        #endregion

        private void GenerateMapMasks(DemoSet set)
        {
            PlaceCoocDataset pcd = new PlaceCoocDataset();

            switch (set)
            {
                case DemoSet.Nunavut:
                    pcd.DrawKML("CanadaEPSG4326.kml", "input/Nunavut_mask", 81.36, 50.38, -139.83, -1.85, 0, 180, 1570, 923, true);
                    break;
                case DemoSet.Bible:
                    pcd.DrawKML("World.kml", "input/MiddleEast_mask", 41.9, 15.35, 12.40, 48.52, 0, 0, 411, 350, true);
                    break;
                case DemoSet.USA:
                    pcd.DrawKML("World.kml", "input/USA_mask", 49.84, 23.08, -126.83, -65.74, 0, 180, 348, 193, true);
                    break;
                default:
                    throw new Exception("Unknown set type");
            }
        }


        private void GenerateHeatmaps(DemoSet set)
        {
            PlaceCoocDataset pcd = new PlaceCoocDataset();

            switch (set)
            {
                case DemoSet.Nunavut:
                    pcd.Init("nunavut.places.txt",
                                "nunavut.coocs.txt", "nunavut.freqs.txt", "Nunavut_mask.png", "output/nunavut/",
                                81.36, 50.38, -139.83, -1.85, 0, 180, 1570, 923);
                    break;

                case DemoSet.Bible:
                    pcd.Init("bible.shortfrequent.chapter.multi.null.places.txt", "bible.shortfrequent.chapter.multi.null.coocs.txt",
                                "bible.shortfrequent.chapter.multi.null.freqs.txt", "MiddleEast_mask.png", "output/bible/",
                                41.9, 15.35, 12.40, 48.52, 0, 0, 411, 350);
                    break;

                case DemoSet.USA:
                    pcd.Init("usa.places.txt", "usa.coocs.txt", "usa.freqs.txt", "USA_mask.png", "output/usa/",
                        49.84, 23.08, -126.83, -65.74, 0, 180, 348, 193);
                    break;

                default:
                    throw new Exception("Unknown set type");
            }
        }


    }
}
