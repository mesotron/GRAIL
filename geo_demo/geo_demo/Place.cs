//FILE: Place.cs

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace geo_demo
{
    class Place
    {
        public string Name { get; set; }
        public string Key { get; set; }
        public List<Coord> Coords = new List<Coord>();
        public bool? Known { get; set; }

        public Place(string name, string key, double lat, double lon, bool? known)
        {
            Name = name;
            Key = key;
            Coords.Add(new Coord(lat, lon));
            Known = known;
        }

        public Place(string name, string key, string coords, bool? known)
        {
            Name = name;
            Key = key;
            string[] toks = coords.Split('|');
            foreach (string t in toks)
            {
                Coords.Add(new Coord(t));
            }
            Known = known;
        }
    }
}
