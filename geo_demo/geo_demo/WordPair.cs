//FILE: WordPair.cs

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Globalization;

public class WordPair
{
    public string Word1 { get; private set; }
    public string Word2 { get; private set; }

    public WordPair(string word1, string word2)
    {
        Word1 = word1;
        Word2 = word2;
    }

    public override string ToString()
    {
        return Word1 + ", " + Word2;
    }

    public override bool Equals(object obj)
    {
        // If parameter is null return false.
        if (obj == null) return false;

        // If parameter cannot be cast to same type, return false.
        WordPair p = obj as WordPair;
        if ((System.Object)p == null) return false;

        // Return true if the fields match:
        return (Word1 == p.Word1) && (Word2 == p.Word2);
    }

    public override int GetHashCode()
    {
        return Word1.GetHashCode() ^ Word2.GetHashCode();
    }
}

