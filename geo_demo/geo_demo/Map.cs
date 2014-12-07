// FILE: Map.cs

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Drawing;

namespace geo_demo
{
    class Map
    {
        private int _width;
        private int _height;
        public int Width { get { return _width; } }
        public int Height { get { return _height; } }

        public List<Point> Points = new List<Point>();
        public Dictionary<Point, Coord> PointsToCoords = new Dictionary<Point, Coord>();
        public bool[,] Area;

        public Map(int width, int height)
        {
            _width = width;
            _height = height;
            Area = new bool[width + 1, height + 1];
        }

        /// <summary>
        /// Constructs a map containing any number of points.
        /// </summary>
        /// <param name="width"></param>
        /// <param name="height"></param>
        /// <param name="point"></param>
        public Map(int width, int height, params Point[] points)
        {
            _width = width;
            _height = height;
            Area = new bool[width + 1, height + 1];
            foreach (Point point in points)
            {
                Points.Add(point);
                Area[point.X, point.Y] = true;
            }
        }

        internal bool IsBoundaryAt(int x, int y)
        {
            if (Area[x, y])
            {
                if (x > 0 && !Area[x - 1, y]) return true;
                if (x < Width - 1 && !Area[x + 1, y]) return true;
                if (y > 0 && !Area[x, y - 1]) return true;
                if (y < Height - 1 && !Area[x, y + 1]) return true;
            }

            return false;
        }
    }
}
