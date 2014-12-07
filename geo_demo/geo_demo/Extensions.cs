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

        //if (word1.CompareTo(word2) > 0)
        //{
        //    Word1 = word2;
        //    Word2 = word1;
        //}
        //else
        //{
            Word1 = word1;
            Word2 = word2;
        //}
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

public static class Extensions
{

	public static Random RNG = new Random();

    public static void WriteSep(this StreamWriter sw, string sep, bool trailingSep, params object[] objs)
    {
        for (int i = 0; i < objs.Length; i++)
        {
            sw.Write(objs[i].ToString());
            if (trailingSep || i != objs.Length - 1) sw.Write(sep);
        }
    }

    public static void WriteSepLine(this StreamWriter sw, string sep, params object[] objs)
    {
        WriteSep(sw, sep, false, objs);
        sw.WriteLine();
    }

    public static void Raise<T>(this EventHandler<T> handler, object sender, T args) where T : EventArgs
    {
        if (handler != null) handler(sender, args);
    }

	public static double? ToNullableDouble(string s)
	{
		double num = 0;

		if (Double.TryParse(s, out num)) {
			return num;
		} else {
			return null;
		}
	}

    public static WordPair Alphabetize(string t1, string t2)
    {
        if (t1.CompareTo(t2) == 1)
            return new WordPair(t2, t1);
        else
            return new WordPair(t1, t2);
    }

    /// <summary>
    /// Returns the minimal element of the given sequence, based on
    /// the given projection.
    /// </summary>
    /// <remarks>
    /// If more than one element has the minimal projected value, the first
    /// one encountered will be returned. This overload uses the default comparer
    /// for the projected type. This operator uses immediate execution, but
    /// only buffers a single result (the current minimal element).
    /// </remarks>
    /// <typeparam name="TSource">Type of the source sequence</typeparam>
    /// <typeparam name="TKey">Type of the projected element</typeparam>
    /// <param name="source">Source sequence</param>
    /// <param name="selector">Selector to use to pick the results to compare</param>
    /// <returns>The minimal element, according to the projection.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="source"/> or <paramref name="selector"/> is null</exception>
    /// <exception cref="InvalidOperationException"><paramref name="source"/> is empty</exception>

    public static TSource MinBy<TSource, TKey>(this IEnumerable<TSource> source,
        Func<TSource, TKey> selector)
    {
        return source.MinBy(selector, Comparer<TKey>.Default);
    }

    /// <summary>
    /// Returns the minimal element of the given sequence, based on
    /// the given projection and the specified comparer for projected values.
    /// </summary>
    /// <remarks>
    /// If more than one element has the minimal projected value, the first
    /// one encountered will be returned. This overload uses the default comparer
    /// for the projected type. This operator uses immediate execution, but
    /// only buffers a single result (the current minimal element).
    /// </remarks>
    /// <typeparam name="TSource">Type of the source sequence</typeparam>
    /// <typeparam name="TKey">Type of the projected element</typeparam>
    /// <param name="source">Source sequence</param>
    /// <param name="selector">Selector to use to pick the results to compare</param>
    /// <param name="comparer">Comparer to use to compare projected values</param>
    /// <returns>The minimal element, according to the projection.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="source"/>, <paramref name="selector"/> 
    /// or <paramref name="comparer"/> is null</exception>
    /// <exception cref="InvalidOperationException"><paramref name="source"/> is empty</exception>

    public static TSource MinBy<TSource, TKey>(this IEnumerable<TSource> source,
        Func<TSource, TKey> selector, IComparer<TKey> comparer)
    {
        if (source == null) throw new ArgumentNullException();
        if (selector == null) throw new ArgumentNullException();
        if (comparer == null) throw new ArgumentNullException();
		
        using (IEnumerator<TSource> sourceIterator = source.GetEnumerator())
        {
            if (!sourceIterator.MoveNext())
            {
                throw new InvalidOperationException("Sequence was empty");
            }
            TSource min = sourceIterator.Current;
            TKey minKey = selector(min);
            while (sourceIterator.MoveNext())
            {
                TSource candidate = sourceIterator.Current;
                TKey candidateProjected = selector(candidate);
                if (comparer.Compare(candidateProjected, minKey) < 0)
                {
                    min = candidate;
                    minKey = candidateProjected;
                }
            }
            return min;
        }
    }

    public static decimal Median<TSource, TKey>(this IEnumerable<TSource> source, Func<TSource, TKey> selector) // IComparer<TKey> comparer)
    {
        if (source == null) throw new ArgumentNullException();
        if (selector == null) throw new ArgumentNullException();
        //if (comparer == null) throw new ArgumentNullException();

        // Create a copy of the input, and sort the copy
        List<TSource> temp = source.ToList();
        temp.Sort();

        int count = temp.Count();
        if (count == 0)
        {
            throw new InvalidOperationException("Empty collection");
        }
        else if (count % 2 == 0)
        {
            // count is even, average two middle elements
            TSource a = temp[count / 2 - 1];
            TSource b = temp[count / 2];
            return (decimal.Parse(a.ToString()) + decimal.Parse(b.ToString())) / 2m;
        }
        else
        {
            // count is odd, return the middle element
            return decimal.Parse(temp[count / 2].ToString());
        }
    }


    public static string ToPercentage(this double d)
    {
        return (d * 100).ToString("N1") + "%";
    }

    public static string[] ParseLineAsCSV(string line, int expectedTokens)
    {
        string[] toReturn = ParseLineAsCSV(line);
        if (toReturn.Length != expectedTokens) SanityCheck.AssertFailed();
        return toReturn;
    }

    public static string[] ParseLineAsCSV(string line)
    {
        List<string> tokens = new List<string>();
        int index = 0;
        bool isEnclosed = false;
        bool justSawQuote = false;
        StringBuilder sb = new StringBuilder();
        while (index < line.Length)
        {
            if (line[index] == ',' && !isEnclosed)
            {
                tokens.Add(sb.ToString());
                sb = new StringBuilder();
                justSawQuote = false;
            }
            else if (line[index] == '"')
            {
                isEnclosed = !isEnclosed;
                if (justSawQuote) sb.Append(line[index]);
            }
            else
            {
                sb.Append(line[index]);
                justSawQuote = false;
            }
        }
        tokens.Add(sb.ToString());
        return tokens.ToArray();
    }

    public static double GetMedian(this IEnumerable<double> source)
    {
        // Create a copy of the input, and sort the copy
        double[] temp = source.ToArray();
        Array.Sort(temp);

        int count = temp.Length;
        if (count == 0)
        {
            throw new InvalidOperationException("Empty collection");
        }
        else if (count % 2 == 0)
        {
            // count is even, average two middle elements
            double a = temp[count / 2 - 1];
            double b = temp[count / 2];
            return (a + b) / 2.0;
        }
        else
        {
            // count is odd, return the middle element
            return temp[count / 2];
        }
    }

    /// <summary>
    /// Compute the distance between two strings.
    /// </summary>
    public static int ComputeLevenshtein(this string s, string t)
    {
        int n = s.Length;
        int m = t.Length;
        int[,] d = new int[n + 1, m + 1];

        // Step 1
        if (n == 0)
        {
            return m;
        }

        if (m == 0)
        {
            return n;
        }

        // Step 2
        for (int i = 0; i <= n; d[i, 0] = i++)
        {
        }

        for (int j = 0; j <= m; d[0, j] = j++)
        {
        }

        // Step 3
        for (int i = 1; i <= n; i++)
        {
            //Step 4
            for (int j = 1; j <= m; j++)
            {
                // Step 5
                int cost = (t[j - 1] == s[i - 1]) ? 0 : 1;

                // Step 6
                d[i, j] = Math.Min(
                    Math.Min(d[i - 1, j] + 1, d[i, j - 1] + 1),
                    d[i - 1, j - 1] + cost);
            }
        }
        // Step 7
        return d[n, m];
    }

    public static int WordCount(this string str)
    {
        return str.Split(new char[] { ' ', '.', '?' },
                         StringSplitOptions.RemoveEmptyEntries).Length;
    }


    /// <summary>
    /// Use the current thread's culture info for conversion
    /// </summary>
    public static string ToTitleCase(this string str)
    {
        var cultureInfo = System.Threading.Thread.CurrentThread.CurrentCulture;
        return cultureInfo.TextInfo.ToTitleCase(str.ToLower());
    }

    public static void Increment<T>(this Dictionary<T, int> source, T key)
    {
        if (source.ContainsKey(key))
        {
            source[key]++;
        }
        else
        {
            source.Add(key, 1);
        }
    }

    // Example uses: see Array<T> below and generalize
    public static IEnumerable<T> IEnum<T>(int length, Func<int, T> initializationFunction)
    {
        return Enumerable.Range(0, length).Select<int, T>(initializationFunction);
    }

    // Example uses:
    //double[] data = Elinq.Array<double>(10, x => rnorm.NextDouble());
    //double[] data = Elinq.Array<double>(10, x => x + 2);
    /*
    public static T[] Array<T>(int length, Func<int, T> initializationFunction)
    {
        return Enumerable.Range(0, length).Select<int, T>(initializationFunction).ToArray<T>();
    }
     */

    // Example uses: see Array<T> above and generalize
    public static List<T> List<T>(int length, Func<int, T> initializationFunction)
    {
        return Enumerable.Range(0, length).Select<int, T>(initializationFunction).ToList<T>();
    }


    /*
     * // Another way to generate sequences; initializes the whole sequence
    //
    // Use: IEnumerable<double> data = Elinq.Init<double>(0d, 100, d => rnorm.NextDouble());
    //
    //
    public static IEnumerable<T> IEnum<T>(T dummyval, int length, Func<T, T> initializationFunction)
    {
        return Enumerable.Repeat(dummyval, length).Select(initializationFunction);
    }
     * 
    // Example use:
    // double[] data = Elinq.Array<double>(0d, 100, f => rnorm.NextDouble());
    public static T[] Array<T>(T dummyval, int length, Func<T, T> initializationFunction)
    {
        return Enumerable.Repeat(dummyval, length).Select(initializationFunction).ToArray<T>();
    }

    // Example use:
    // List<double> data = Elinq.List<double>(0d, 100, f => rnorm.NextDouble());
    public static List<T> List<T>(T dummyval, int length, Func<T, T> initializationFunction)
    {
        return Enumerable.Repeat(dummyval, length).Select(initializationFunction).ToList<T>();
    }
     */

    //
    // Generating sequences
    //
    // Use: IEnumerable<double> data = Elinq.Generate<double>(() => (double?)rnorm.NextDouble());
    //
    public static IEnumerable<T> Generate<T>(Func<T> generator) where T : class
    {
        if (generator == null) throw new ArgumentNullException("generator");

        T t;
        while ((t = generator()) != null)
        {
            yield return t;
        }
    }

    public static IEnumerable<T> Generate<T>(Func<Nullable<T>> generator) where T : struct
    {
        if (generator == null) throw new ArgumentNullException("generator");

        Nullable<T> t;
        while ((t = generator()).HasValue)
        {
            yield return t.Value;
        }
    }

    public static IEnumerable<T> FromEnumerator<T>(IEnumerator<T> enumerator)
    {
        if (enumerator == null) throw new ArgumentNullException("enumerator");

        while (enumerator.MoveNext())
        {
            yield return enumerator.Current;
        }
    }

    public static IEnumerable<T> Single<T>(T value)
    {
        return Enumerable.Repeat(value, 1);
    }

    //
    // I/O
    //
    public static IEnumerable<string> ReadLinesFromFile(string path)
    {
        if (path == null) throw new ArgumentNullException("path");
        using (StreamReader file = new StreamReader(path))
        {
            string line;
            while ((line = file.ReadLine()) != null) yield return line;
        }
    }

    public static IEnumerable<string> ReadLinesFromConsole()
    {
        return ReadLinesFrom(Console.In);
    }

    public static IEnumerable<string> ReadLinesFrom(TextReader reader)
    {
        if (reader == null) throw new ArgumentNullException("reader");

        return Generate(() => reader.ReadLine());
    }

    public static void WriteLinesTo<T>(this IEnumerable<T> lines, TextWriter writer)
    {
        if (lines == null) throw new ArgumentNullException("lines");
        if (writer == null) throw new ArgumentNullException("writer");

        lines.ForEach((line) => writer.WriteLine(line.ToString()));
    }

    public static void WriteLinesToConsole<T>(this IEnumerable<T> lines)
    {
        lines.WriteLinesTo(Console.Out);
    }

    public static void WriteLinesToFile<T>(this IEnumerable<T> lines, string path)
    {
        if (path == null) throw new ArgumentNullException("path");

        using (TextWriter file = new StreamWriter(path))
        {
            lines.WriteLinesTo(file);
        }
    }

    //
    // Side effects
    //
    public static IEnumerable<T> Do<T>(this IEnumerable<T> source, Action<T> action)
    {
        if (source == null) throw new ArgumentNullException("source");
        if (action == null) throw new ArgumentNullException("action");

        foreach (T elem in source)
        {
            action(elem);
            yield return elem;
        }
    }

    public static void ForEach<T>(this IEnumerable<T> source, Action<T> action)
    {
        if (source == null) throw new ArgumentNullException("source");
        if (action == null) throw new ArgumentNullException("action");

        foreach (T elem in source)
        {
            action(elem);
        }
    }

    //
    // ToStringPretty
    //
    public static string ToStringPretty<T>(this IEnumerable<T> source)
    {
        return ToStringPretty(source, ",");
    }

    public static string ToStringPretty<T>(this IEnumerable<T> source, string delimiter)
    {
        return ToStringPretty(source, "", delimiter, "");
    }

    public static string ToStringPretty<T>(this IEnumerable<T> source, string before, string delimiter, string after)
    {
        if (source == null) throw new ArgumentNullException("source");

        StringBuilder result = new StringBuilder();
        result.Append(before);

        bool firstElement = true;
        foreach (T elem in source)
        {
            if (firstElement) firstElement = false;
            else result.Append(delimiter);

            result.Append(elem.ToString());
        }

        result.Append(after);
        return result.ToString();
    }

    //
    // Other
    //
    public static IEnumerable<TOut> Combine<TIn1, TIn2, TOut>(
        this IEnumerable<TIn1> in1, IEnumerable<TIn2> in2, Func<TIn1, TIn2, TOut> func)
    {
        if (in1 == null) throw new ArgumentNullException("in1");
        if (in2 == null) throw new ArgumentNullException("in2");
        if (func == null) throw new ArgumentNullException("func");

        using (var e1 = in1.GetEnumerator())
        using (var e2 = in2.GetEnumerator())
        {
            while (e1.MoveNext() && e2.MoveNext())
            {
                yield return func(e1.Current, e2.Current);
            }
        }
    }


    public static IEnumerable<T> Shuffle<T>(this IEnumerable<T> source)
    {
        if (source == null) throw new ArgumentNullException("source");

        T[] array = source.ToArray();

        for (int i = 0; i < array.Length; i++)
        {
            int r = RNG.Next(i + 1);
            T tmp = array[r];
            array[r] = array[i];
            array[i] = tmp;
        }

        return array;
    }
}
