using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Windows.Forms;

namespace geo_demo
{
    public partial class PlaceForm : Form
    {
        public PlaceForm()
        {
            InitializeComponent();
        }

        private void PlaceForm_Load(object sender, EventArgs e)
        {
            DataTable table1 = new DataTable("People");
            table1.Columns.Add("id");
            table1.Columns.Add("Name");
            table1.Columns.Add("Age");
            table1.Rows.Add(1, "Jack", 18);
            table1.Rows.Add(2, "Tim", 18);

            DataSet set = new DataSet("SetPeople");
            set.Tables.Add(table1);

            dataGridView1.DataSource = set;
            dataGridView1.DataSource = table1;

            dataGridView1.Update();

            
        }
    }
}
