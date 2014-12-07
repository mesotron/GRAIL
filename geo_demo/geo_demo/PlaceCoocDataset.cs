//FILE: PlaceCoocDataset.cs

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Text.RegularExpressions;
using System.Drawing;
using System.Drawing.Drawing2D;
using System.Drawing.Imaging;

namespace geo_demo
{
    class PlaceCoocDataset
    {
        const bool drawMaps = true;

        SortedDictionary<string, Place> _places = new SortedDictionary<string, Place>();
        Dictionary<WordPair, int> _coocs = new Dictionary<WordPair, int>();      // <known, unknown>
        Dictionary<WordPair, double> _pmis = new Dictionary<WordPair, double>();      // <known, unknown>
        SortedDictionary<string, int> _freqs = new SortedDictionary<string, int>();
        HashSet<string> _knownSet = new HashSet<string>();
        HashSet<string> _unknownSet = new HashSet<string>();

        SortedDictionary<string, Map> _maps = new SortedDictionary<string, Map>();
        Dictionary<string, Dictionary<Point, double>> _toponymsToLLs = new Dictionary<string, Dictionary<Point, double>>();

        public double TopLat { get; set; }
        public double BottomLat { get; set; }
        public double LeftLon { get; set; }
        public double RightLon { get; set; }

        public double LatOffset { get; set; }
        public double LonOffset { get; set; }

        public int MapWidth { get; set; }
        public int MapHeight { get; set; }



        /// <summary>
        /// Load a complete dataset of places and co-occurrences.
        /// 
        /// Acceptable placesFile format:
        /// Must have this header: name\tkey\tlat\tlon\tknown
        /// 
        /// Acceptable co-occurrence grid formats:
        /// - unknowns in cols, knowns in rows
        /// OR
        /// - everything is null ("cross-validated") and the rows are the same as the cols
        /// Either way, both rows and cols should be labeled
        /// 
        /// Acceptable freqFile format:
        /// key\tfreq
        /// 
        /// </summary>
        public void Init(string placesFile, string coocsFile, string freqFile, string mainlandBmpFile, string outputPath,
            double topLat, double bottomLat, double leftLon, double rightLon,
            double latOffset, double lonOffset, int outputWidth, int outputHeight)
        {
            string inputPath = "input/";
            string outputBestFile = Path.GetFileNameWithoutExtension(placesFile) + ".best.txt";

            topLat += latOffset;
            bottomLat += latOffset;
            leftLon += lonOffset;
            rightLon += lonOffset;

            TopLat = topLat;
            BottomLat = bottomLat;
            LeftLon = leftLon;
            RightLon = rightLon;
            LatOffset = latOffset;
            LonOffset = lonOffset;
            MapWidth = outputWidth;
            MapHeight = outputHeight;

            // Load places, co-occurrences, and freqs
            LoadPlaces(inputPath + placesFile);
            LoadCoocs(inputPath + coocsFile);
            LoadFreqs(inputPath + freqFile);

            // Start the computation proper
            ComputePMIs(_coocs, _freqs);


            Dictionary<string, Point> unknownsToPoints;     // Maps each unknown place name to the point on the map that corresponds to its actual lat/lon.

            Maths.MercatorInfo mi = InitializeMaps(out unknownsToPoints, topLat, bottomLat, leftLon, rightLon, latOffset, lonOffset, outputWidth, outputHeight);
            Map mainland = ReadMap(new Bitmap(inputPath + mainlandBmpFile), Color.White, mi);


            // Continue computation as in geo_confidence.LLMapper.Init()
            // Compute shortest distances from points on raster map to known/training villages
            Dictionary<Point, List<float>> pointsToDistanceVectors = ComputeShortestDistances(mainland);

            // Compute probability of interaction--according to gravity model--between {an imaginary village at each point} and {each actual village} at beta of .2
            Dictionary<Point, List<float>> pointsToExpectedProbabilityVectors = ComputeProbabilityVectors(pointsToDistanceVectors, 0.2);

            // 1a) For each village, for each point: calculate the log-likelihood of observing the data that were actually observed
            // 1b) Output the "best" lat-lon coordinate, along with where the place actually is
            //
            _toponymsToLLs = ComputeLLs(outputPath, outputBestFile, mainland, pointsToDistanceVectors, pointsToExpectedProbabilityVectors);

            // 2) Figure out the appropriate mapping between log-likelihoods and colors, based on observed probabilities
            SortedDictionary<int, double> percentsToNormalizedLLThresholds = ComputeNormalizedLLThresholds(ref unknownsToPoints);
            Dictionary<int, Color> percentsToColors = GetColorMapping();
            SaveColorKey(outputPath, percentsToColors, percentsToNormalizedLLThresholds);

            // 3) Figure out how much time we saved a hypothetical analyst
            ComputeTimeSaved(outputPath, unknownsToPoints);

            // 4) Make the maxes and mins the ends of a color spectrum, and draw log-likelihood maps accordingly

            if (drawMaps)
            {
                foreach (string toponym in _unknownSet)
                {
                    var lls = _toponymsToLLs[toponym];

                    double min = lls.Where(x => !(double.IsNaN(x.Value) || double.IsInfinity(x.Value))).Min(x => x.Value);
                    double max = lls.Max(x => x.Value);
                    double range = max - min;

                    Bitmap b = new Bitmap(MapWidth, MapHeight);
                    BitmapData bmData = b.LockBits(new Rectangle(0, 0, b.Width, b.Height),
                            ImageLockMode.ReadWrite, PixelFormat.Format24bppRgb);
                    int stride = bmData.Stride;
                    System.IntPtr Scan0 = bmData.Scan0;
                    unsafe
                    {
                        byte* p = (byte*)(void*)Scan0;
                        int nOffset = stride - b.Width * 3;

                        for (int y = 0; y < b.Height; ++y)
                        {
                            for (int x = 0; x < b.Width; ++x)
                            {
                                if (!mainland.Area[x, y])
                                {
                                    p[0] = 0;
                                    p[1] = 0;
                                    p[2] = 0;
                                }
                                else
                                {
                                    double ll = lls[new Point(x, y)];

                                    double normalized = (ll - min) / range;

                                    Color c = NormalizedLLToColor(normalized, percentsToNormalizedLLThresholds, percentsToColors);
                                    p[0] = c.B;
                                    p[1] = c.G;
                                    p[2] = c.R;
                                }

                                p += 3;

                            }
                            p += nOffset;
                        }
                    }
                    b.UnlockBits(bmData);
                    b.Save(outputPath + toponym + ".png", ImageFormat.Png);

                    //Draw the actual location of the place if it is known

                    Graphics g = Graphics.FromImage(b);
                    Point pt = unknownsToPoints[toponym];

                    g.FillRectangle(Brushes.LimeGreen, new Rectangle(pt.X - 2, pt.Y - 2, 4, 4));
                    g.FillRectangle(Brushes.Black, new Rectangle(pt.X - 1, pt.Y - 1, 2, 2));


                    b.Save(outputPath + toponym + "_answer.png", ImageFormat.Png);

                }
            }

           

        }



        public void DrawKML(string kmlFile, string partOfWorld, double topLat, double bottomLat, double leftLon, double rightLon, double latOffset, double lonOffset, int outputWidth, int outputHeight, bool setNegValuesToZero)
        {
            topLat += latOffset;
            bottomLat += latOffset;
            leftLon += lonOffset;
            rightLon += lonOffset;

            // Load KML file

            var multigeometries = new List<List<List<Coord>>>();
            using (StreamReader sr = new StreamReader(kmlFile))
            {
                while (!sr.EndOfStream)
                {
                    string line = sr.ReadLine().Trim();
                    if (line == "") continue;

                    if (line.StartsWith("<MultiGeometry>") || line.StartsWith("<Polygon>"))
                    {
                        List<List<Coord>> polygons = new List<List<Coord>>();

                        line = line.Replace("<Polygon>", "|");
                        line = line.Replace("</Polygon>", "|");
                        line = line.Replace("<outerBoundaryIs><LinearRing><coordinates>", "|");
                        line = line.Replace("</coordinates></LinearRing></outerBoundaryIs>", "|");
                        line = line.Replace("<innerBoundaryIs><LinearRing><coordinates>", "|");
                        line = line.Replace("</coordinates></LinearRing></innerBoundaryIs>", "|");
                        line = line.Replace("||", "|");
                        line = line.Replace("<MultiGeometry>|", "");
                        line = line.Replace("|</MultiGeometry>", "");

                        string[] polygonLines = line.Split('|');

                        foreach (string p in polygonLines)
                        {
                            List<Coord> polygon = new List<Coord>();
                            string[] coords = p.Split(' ');
                            foreach (string c in coords)
                            {
                                if (c != "")
                                {
                                    string[] toks = c.Split(',');
                                    double lon = double.Parse(toks[0]) + lonOffset; // lon comes first in kml
                                    double lat = double.Parse(toks[1]) + latOffset;

                                    if (setNegValuesToZero)
                                    {
                                        if (lon < 0) lon = 0;
                                        if (lat < 0) lat = 0;
                                    }
                                    else
                                    {
                                        if (lon < 0) SanityCheck.AssertFailed();    // We have negative longitudes: try a different lonOffset
                                        if (lat < 0) SanityCheck.AssertFailed();    // We have negative latitudes: try a different latOffset
                                    }

                                    polygon.Add(new Coord(lat, lon));
                                }
                            }
                            polygons.Add(polygon);
                        }
                        multigeometries.Add(polygons);
                    }
                }
            }

            // Create bitmap & graphics device to draw to
            Bitmap b = new Bitmap(outputWidth, outputHeight);
            Graphics g = Graphics.FromImage(b);
            g.Clear(Color.Black);


            Maths.MercatorInfo mi = new Maths.MercatorInfo(bottomLat.ToRadians(), topLat.ToRadians(), leftLon.ToRadians(), rightLon.ToRadians(), outputWidth, outputHeight);

            foreach (List<List<Coord>> polygons in multigeometries)
            {
                foreach (List<Coord> polygon in polygons)
                {
                    if (polygon.Count() > 0)
                    {
                        PointF[] points = polygon.Select(x => Maths.ProjectToMercator(x, mi)).ToArray();
                        g.FillPolygon(Brushes.White, points);

                        for (int i = 0; i < points.Length; i++)
                        {
                            Coord projectedBack = Maths.ProjectBackFromMercator(points[i], mi);

                        }
                    }
                }
            }

            b.Save(partOfWorld + ".png", ImageFormat.Png);
        }



        private void ComputeTimeSaved(string outputPath, Dictionary<string, Point> unknownsToPoints)
        {
            
            using (StreamWriter sw = new StreamWriter(outputPath + "time_saved.txt"))
            {
                sw.WriteSepLine("\t", "key", "time_if_systematic", "time_if_random", "time_saved");

                foreach (string key in unknownsToPoints.Keys)
                {
                    Point pt = unknownsToPoints[key];
                    double ll = _toponymsToLLs[key][pt];

                    int totalPoints = _toponymsToLLs[key].Count();
                    int timeIfRandom = totalPoints / 2;
                    int timeIfSystematic = _toponymsToLLs[key].Where(x => x.Value >= ll).Count();
                    double percentTimeSaved = 1.0 - (timeIfSystematic / (double)timeIfRandom);

                    sw.WriteSepLine("\t", key, timeIfSystematic, timeIfRandom, percentTimeSaved.ToPercentage());
                }
            }
            
        }


        private void SaveColorKey(string outputPath, Dictionary<int, Color> percentsToColors, SortedDictionary<int, double> percentsToNormalizedLLThresholds)
        {
            // Save color key (mapping from percents to colors)
            //
            Font font = new Font("Arial", 8, FontStyle.Regular);
            Bitmap colorkey = new Bitmap(2000, 30);
            Graphics g = Graphics.FromImage(colorkey);

            for (int i = 1; i < 100; i++)
            {
                g.FillRectangle(new SolidBrush(percentsToColors[i]), new Rectangle(i * 20, 0, 20, 10));
                g.DrawString(i.ToString(), font, Brushes.Black, i * 20, 15);
            }
            colorkey.Save(outputPath + "colorkey.png");

            // Save mapping from percents to normalized ll thresholds
            //
            using (StreamWriter sw = new StreamWriter(outputPath + "percentsToNormalizedLLThresholds.txt"))
            {
                for (int i = 1; i < 100; i++)
                {
                    sw.WriteSepLine("\t", i, percentsToNormalizedLLThresholds[i].ToString("N5"));
                }
            }
        }

        private SortedDictionary<int, double> ComputeNormalizedLLThresholds(ref Dictionary<string, Point> unknownsToPoints)
        {
            SortedDictionary<int, double> percentsToNormalizedLLThresholds = new SortedDictionary<int, double>();
            // ^^^ How to interpret: "KEY% of villages have normalized log likelihoods higher than or equal to VALUE"

            // Get LL at the actual location of each unknown/cross-validated village
            List<double> normalizedLLs = new List<double>();
            foreach (string unknown in _unknownSet)
            {
                Point pt = unknownsToPoints[unknown];
                double ll;

                if (_toponymsToLLs[unknown].ContainsKey(pt))
                {
                    ll = _toponymsToLLs[unknown][pt];
                }
                else
                {
                    // This can happen if the actual location is slightly in the sea, or slightly outside of the window
                    // In this case, find the closest land point that we do have a LL for
                    double minDistance = double.MaxValue;
                    Point bestPoint = new Point(0, 0);
                    foreach (Point pt2 in _toponymsToLLs[unknown].Keys)
                    {
                        double dist = Maths.EucDistance(pt.X, pt.Y, pt2.X, pt2.Y);
                        if (dist < minDistance)
                        {
                            minDistance = dist;
                            bestPoint = pt2;
                        }
                    }

                    unknownsToPoints[unknown] = bestPoint;
                    ll = _toponymsToLLs[unknown][bestPoint];
                }

                double min = _toponymsToLLs[unknown].Where(x => !(double.IsNaN(x.Value) || double.IsInfinity(x.Value))).Min(x => x.Value);
                double max = _toponymsToLLs[unknown].Max(x => x.Value);
                double range = max - min;

                double normalized = (ll - min) / range;
                normalizedLLs.Add(normalized);
            }

            // Sort LLs from largest to smallest, and iterate through to find the HIGHEST (smallest area) LL thresholds t-sub-n such that
            // * 1% of actual locations have LLs >= t1
            // * 2% of actual locations have LLs >= t2
            // etc.

            normalizedLLs.Sort();
            normalizedLLs.Reverse();

            double interval = 1.0 / normalizedLLs.Count();

            for (int i = 0; i < normalizedLLs.Count(); i++)
            {
                double ll = normalizedLLs[i];
                int min = (int)Math.Round(interval * i * 100);
                int max = (int)Math.Round(interval * (i + 1) * 100);
                for (int j = min; j <= max; j++)
                {
                    if (!percentsToNormalizedLLThresholds.ContainsKey(j))
                    {
                        percentsToNormalizedLLThresholds.Add(j, ll);
                    }
                }
            }

            return percentsToNormalizedLLThresholds;
        }

        private Dictionary<int, Color> GetColorMapping()
        {
            Dictionary<int, Color> colorMapping = new Dictionary<int, Color>();

            float interval = 240 / 100f;

            for (int i = 0; i <= 100; i++)
			{
                Color c = ColorFromAhsb(255, interval * i, 1f, .5f);
                colorMapping.Add(i, c);
			}

            return colorMapping;
        }

        private Color NormalizedLLToColor(double ll, SortedDictionary<int, double> percentsToLLThresholds, Dictionary<int, Color> percentsToColors)
        {
            int percent = -1;
            foreach (KeyValuePair<int, double> kv in percentsToLLThresholds)
            {
                if (ll >= kv.Value)
                {
                    percent = kv.Key;
                    break;
                }
            }
            if (percent == -1) percent = 100;
            return percentsToColors[percent];
        }


        // 1a) For each village, for each point: calculate the log-likelihood of observing the data that were actually observed
        // 1b) Output the "best" lat-lon coordinate, along with where the place actually is
        //
        private Dictionary<string, Dictionary<Point, double>> ComputeLLs(string outputPath, string outputBestFile, Map mainland, Dictionary<Point, List<float>> pointsToDistanceVectors, Dictionary<Point, List<float>> pointsToExpectedProbabilityVectors)
        {
            pointsToDistanceVectors = null;
            GC.Collect();

            var knownList = _maps.Keys.ToList();

            using (StreamWriter sw = new StreamWriter(outputPath + outputBestFile))
            {
                sw.WriteLine("id\tname\tlat\tlon\testlat\testlon");

                int i = 0;
                foreach (string unknown in _unknownSet)
                {
                    int indexToOmit = -1;
                    if (knownList.Contains(unknown)) indexToOmit = knownList.IndexOf(unknown);

                    _toponymsToLLs.Add(unknown, new Dictionary<Point, double>());
                    double maxLL = Double.MinValue;
                    Coord bestCoord = new Coord(0, 0);

                    for (int x = 0; x < MapWidth; x++)
                    {
                        for (int y = 0; y < MapHeight; y++)
                        {
                            if (mainland.Area[x, y])
                            {
                                Point pt1 = new Point(x, y);
                                double ll = 0;
                                List<float> normalizedVector = Renormalize(pointsToExpectedProbabilityVectors[pt1], indexToOmit);

                                int j = 0;
                                foreach (string known in _maps.Keys)
                                {
                                    if (unknown != known) // We are not allowed any knowledge about the actual location of the place.
                                    {
                                        ll += _pmis[new WordPair(known, unknown)] * Math.Log(normalizedVector[j]);
                                    }
                                    j++;
                                }

                                // We now have the LL for "unknown", for this point; store it
                                if (drawMaps) _toponymsToLLs[unknown].Add(pt1, ll);

                                // If it is better than the best LL so far obtained, save it
                                if (ll > maxLL)
                                {
                                    maxLL = ll;
                                    bestCoord = mainland.PointsToCoords[pt1];
                                }
                            }
                        }
                    }

                    // Write the best coordinate obtained for this unknown
                    bestCoord.Lat -= LatOffset;
                    bestCoord.Lon -= LonOffset;

                    if (_places[unknown].Coords.Count == 1)
                    {
                        Coord actual = _places[unknown].Coords[0];   // We can do this in this case since we know each place is represented by only one point
                        sw.WriteLine(i.ToString() + "\t" + unknown + "\t" + actual.Lat.ToString() + "\t" + actual.Lon.ToString() + "\t" + bestCoord.Lat.ToString() + "\t" + bestCoord.Lon.ToString());
                    }
                    else
                    {
                        sw.WriteLine(i.ToString() + "\t" + unknown + "\t" + String.Join("|", _places[unknown].Coords) + "\t" + bestCoord.Lat.ToString() + "\t" + bestCoord.Lon.ToString());
                    }
                    sw.Flush();
                    i++;
                }

            }
            return _toponymsToLLs;
        }


        private List<float> Renormalize(List<float> vector, int indexToOmit)
        {
            if (indexToOmit == -1)
            {
                return vector;
            }
            else
            {
                var vector2 = new List<float>(vector);
                vector2[indexToOmit] = 0;
                double oneOverSum = 1.0 / vector2.Sum();
                return vector2.Select(x => (float)(x * oneOverSum)).ToList();
            }
        }

        private Dictionary<Point, List<float>> ComputeProbabilityVectors(Dictionary<Point, List<float>> pointsToDistanceVectors, double beta)
        {
            var pointsToProbabilityVectors = new Dictionary<Point, List<float>>();

            // Calculate probability vectors
            foreach (Point pt in pointsToDistanceVectors.Keys)
            {
                var probabilityVector = new List<float>();
                for (int i = 0; i < _maps.Keys.Count(); i++)
                {
                    //This is what you would do if you wanted to normalize by population
                    //double weight = abbrsToPopulations[abbr] / Math.Pow(pointsToDistanceVectors[pt][i], beta);

                    //Since we are already essentially normalizing for frequency by using PMI, treat populations as equivalent
                    float weight = (float)(1.0 / Math.Pow(pointsToDistanceVectors[pt][i], beta));

                    probabilityVector.Add(weight);
                }

                double sum = probabilityVector.Sum();
                pointsToProbabilityVectors.Add(pt, probabilityVector.Select(x => (float)(x / sum)).ToList());
            }

            return pointsToProbabilityVectors;
        }



        private Dictionary<Point, List<float>> ComputeShortestDistances(Map mainland)
        {
            var pointsToDistanceVectors = new Dictionary<Point, List<float>>();

            for (int x = 0; x < MapWidth; x++)
            {
                for (int y = 0; y < MapHeight; y++)
                {
                    if (mainland.Area[x, y])
                    {
                        Point pt1 = new Point(x, y);
                        List<float> distList = new List<float>();
                        foreach (string key in _maps.Keys)
                        {
                            distList.Add((float)GetShortestGISDistance(mainland, pt1, key));
                        }
                        pointsToDistanceVectors.Add(pt1, distList);

                    }
                }
            }

            return pointsToDistanceVectors;
        }

        /// <summary>
        /// Get shortest GIS distance from a point to a place (using coordinates, before projection).
        /// </summary>
        /// <param name="pt1"></param>
        /// <param name="knownPlace"></param>
        /// <returns></returns>
        private double GetShortestGISDistance(Map mainland, Point pt1, string knownPlace)
        {
            List<double> distances = new List<double>();
            Coord c1 = mainland.PointsToCoords[pt1];

            foreach (Point pt2 in _maps[knownPlace].Points)
            {
                Coord c2 = _maps[knownPlace].PointsToCoords[pt2];
                if (c1.Equals(c2)) return 0;
                distances.Add(Maths.GISDistance(c1.Lat, c1.Lon, c2.Lat, c2.Lon));
            }

            return distances.Min();
        }

        private Maths.MercatorInfo InitializeMaps(out Dictionary<string, Point> unknownsToPoints, double topLat, double bottomLat, double leftLon, double rightLon, double latOffset, double lonOffset, int outputWidth, int outputHeight)
        {
            // Initialize maps directly from latitudes and longitudes already loaded in
            // We are initializing maps for KNOWN places and TRAINING (cross-validated) places.
            //
            Maths.MercatorInfo mi = new Maths.MercatorInfo(bottomLat.ToRadians(), topLat.ToRadians(), leftLon.ToRadians(), rightLon.ToRadians(), outputWidth, outputHeight);
            unknownsToPoints = new Dictionary<string, Point>();

            foreach (KeyValuePair<string, Place> kv in _places)
            {
                Place p = kv.Value;

                foreach (Coord c in p.Coords)
                {
                    Coord coord = new Coord(c.Lat + latOffset, c.Lon + lonOffset);
                    PointF ptF = Maths.ProjectToMercator(coord, mi);
                    Point pt = new Point((int)Math.Round(ptF.X), (int)Math.Round(ptF.Y));

                    if (p.Known == false || p.Known == null)
                    {
                        unknownsToPoints.Add(p.Key, pt);
                    }

                    if (p.Known == true || p.Known == null)
                    {
                        Map m = new Map(MapWidth, MapHeight);
                        m.Points.Add(pt);
                        m.PointsToCoords.Add(pt, coord);
                        _maps.Add(kv.Key, m);
                    }
                }
            }

            if (_maps.Count() == 0) SanityCheck.AssertFailed();
            return mi;
        }



        Map ReadMap(Bitmap b, Color color, Maths.MercatorInfo mi)
        {
            Map map = new Map(MapWidth, MapHeight);

            BitmapData bmData = b.LockBits(new Rectangle(0, 0, b.Width, b.Height),
                ImageLockMode.ReadWrite, PixelFormat.Format24bppRgb);
            int stride = bmData.Stride;
            System.IntPtr Scan0 = bmData.Scan0;
            unsafe
            {
                byte* p = (byte*)(void*)Scan0;
                int nOffset = stride - b.Width * 3;
                byte red, green, blue;

                for (int y = 0; y < b.Height; ++y)
                {
                    for (int x = 0; x < b.Width; ++x)
                    {
                        blue = p[0];
                        green = p[1];
                        red = p[2];

                        if (color.R == red && color.G == green && color.B == blue)
                        {
                            Point pt = new Point(x, y);
                            Coord c = Maths.ProjectBackFromMercator(new PointF(x, y), mi);

                            map.Points.Add(pt);
                            map.PointsToCoords.Add(pt, c);
                            map.Area[x, y] = true;
                        }
                        else
                        {
                            p[0] = 0;
                            p[1] = 0;
                            p[2] = 0;
                        }

                        p += 3;

                    }
                    p += nOffset;
                }
            }
            b.UnlockBits(bmData);

            return map;
        }


        private void ComputePMIs(Dictionary<WordPair, int> coocs, SortedDictionary<string, int> freqs)
        {
            // Convert coocs to PMIs
            //
            foreach (string k in _knownSet)
            {
                foreach (string u in _unknownSet)
                {
                    WordPair wp = new WordPair(k, u);

                    if (!_pmis.ContainsKey(wp))
                    {
                        double denom = freqs[k] * (freqs[u] * .0001);
                        if (denom < 0) SanityCheck.AssertFailed();
                        _pmis.Add(wp, coocs[wp] / denom);
                    }
                }
            }
        }

        private void LoadFreqs(string freqFile)
        {
            // Load frequencies
            //
            using (StreamReader sr = new StreamReader(freqFile))
            {
                string header = sr.ReadLine().Trim();
                if (header != "key\tfreq") SanityCheck.AssertFailed();

                while (!sr.EndOfStream)
                {
                    string line = sr.ReadLine().Trim();
                    if (line == "") continue;
                    string[] tokens = line.Split('\t');

                    if (!_places.ContainsKey(tokens[0])) SanityCheck.AssertFailed();
                    _freqs.Add(tokens[0], int.Parse(tokens[1]));
                }
            }
        }

        private void LoadCoocs(string coocsFile)
        {
            // Load coocs
            //
            using (StreamReader sr = new StreamReader(coocsFile))
            {
                string header = sr.ReadLine().Trim();
                string[] cols = header.Split('\t');

                foreach (string col in cols)
                {
                    if (!_places.ContainsKey(col)) SanityCheck.AssertFailed();
                    if (_places[col].Known == true) SanityCheck.AssertFailed();  // Cols should be unknown or null
                }

                while (!sr.EndOfStream)
                {
                    string line = sr.ReadLine().Trim();
                    if (line == "") continue;
                    string[] tokens = line.Split('\t');

                    string rowLabel = tokens[0];
                    if (!_places.ContainsKey(rowLabel)) SanityCheck.AssertFailed();
                    if (_places[rowLabel].Known == false) SanityCheck.AssertFailed();    // Rows should be known or null

                    for (int i = 1; i < tokens.Count(); i++)
                    {
                        WordPair wp = new WordPair(rowLabel, cols[i - 1]);
                        if (!_coocs.ContainsKey(wp))
                        {
                            _coocs.Add(wp, int.Parse(tokens[i]));
                        }
                    }
                }
            }
        }

        private void LoadPlaces(string placesFile)
        {
            // Load places
            //
            using (StreamReader sr = new StreamReader(placesFile))
            {
                string header = sr.ReadLine().Trim();
                if (header != "name\tkey\tlat\tlon\tknown" && header != "name\tkey\tcoords\tknown") SanityCheck.AssertFailed();

                while (!sr.EndOfStream)
                {
                    string line = sr.ReadLine().Trim();
                    if (line == "") continue;
                    string[] tokens = line.Split('\t');

                    bool? known = null;
                    string knownTok = tokens[tokens.Length - 1].ToLower();
                    if (knownTok != "cross-validated")
                    {
                        if (knownTok == "true" || knownTok == "known")
                        {
                            known = true;
                        }
                        else if (knownTok == "false" || knownTok == "unknown")
                        {
                            known = false;
                        }
                        else
                        {
                            SanityCheck.AssertFailed();
                        }
                    }

                    if (tokens.Length == 5)
                    {
                        Place p = new Place(tokens[0], tokens[1], double.Parse(tokens[2]), double.Parse(tokens[3]), known);
                        _places.Add(p.Key, p);
                    }
                    else if (tokens.Length == 4)
                    {
                        Place p = new Place(tokens[0], tokens[1], tokens[2], known);
                        _places.Add(p.Key, p);
                    }
                }
            }

            // Sanity check to make sure that places are all bools (true/false) and that there is at least one unknown place, or else all null
            int nullCount = 0;
            int unkCount = 0;
            foreach (Place p in _places.Values)
            {
                if (p.Known == null)
                    nullCount++;
                else if (p.Known == false)
                    unkCount++;
            }
            if (nullCount == 0)
            {
                if (unkCount == 0) SanityCheck.AssertFailed();
            }
            else
            {
                if (nullCount != _places.Count()) SanityCheck.AssertFailed();
            }

            // Put copies of keys in the "knowns" and "unknowns" lists
            foreach (Place p in _places.Values)
            {
                if (p.Known == null || p.Known == true)
                {
                    _knownSet.Add(p.Key);
                }
                if (p.Known == null || p.Known == false)
                {
                    _unknownSet.Add(p.Key);
                }
            }
        }

        public static Color ColorFromAhsb(int a, float h, float s, float b)
        {

            if (0 > a || 255 < a)
            {
                SanityCheck.AssertFailed();
            }
            if (0f > h || 360f < h)
            {
                SanityCheck.AssertFailed();
            }
            if (0f > s || 1f < s)
            {
                SanityCheck.AssertFailed();
            }
            if (0f > b || 1f < b)
            {
                SanityCheck.AssertFailed();
            }

            if (0 == s)
            {
                return Color.FromArgb(a, Convert.ToInt32(b * 255),
                  Convert.ToInt32(b * 255), Convert.ToInt32(b * 255));
            }

            float fMax, fMid, fMin;
            int iSextant, iMax, iMid, iMin;

            if (0.5 < b)
            {
                fMax = b - (b * s) + s;
                fMin = b + (b * s) - s;
            }
            else
            {
                fMax = b + (b * s);
                fMin = b - (b * s);
            }

            iSextant = (int)Math.Floor(h / 60f);
            if (300f <= h)
            {
                h -= 360f;
            }
            h /= 60f;
            h -= 2f * (float)Math.Floor(((iSextant + 1f) % 6f) / 2f);
            if (0 == iSextant % 2)
            {
                fMid = h * (fMax - fMin) + fMin;
            }
            else
            {
                fMid = fMin - h * (fMax - fMin);
            }

            iMax = Convert.ToInt32(fMax * 255);
            iMid = Convert.ToInt32(fMid * 255);
            iMin = Convert.ToInt32(fMin * 255);

            switch (iSextant)
            {
                case 1:
                    return Color.FromArgb(a, iMid, iMax, iMin);
                case 2:
                    return Color.FromArgb(a, iMin, iMax, iMid);
                case 3:
                    return Color.FromArgb(a, iMin, iMid, iMax);
                case 4:
                    return Color.FromArgb(a, iMid, iMin, iMax);
                case 5:
                    return Color.FromArgb(a, iMax, iMin, iMid);
                default:
                    return Color.FromArgb(a, iMax, iMid, iMin);
            }
        }

    }




}
