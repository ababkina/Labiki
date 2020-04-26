using System;

namespace lab4
{
    public class Matrix
    { ////////////////
        public static double[,] MatrixDecompose(double[,] matrix, out int[] perm, out int toggle)
        {
            // Разложение LUP Дулитла. Предполагается,
            // что матрица квадратная.
            int n = 6; // для удобства
            double[][] result = new double[n][];
            for (int i = 0; i < n; i++)
            {
                result[i] = new double[n];
                for (int j = 0; j < n; j++)
                    result[i][j] = matrix[i, j];
            }
            perm = new int[n];
            for (int i = 0; i < n; ++i) { perm[i] = i; }
            toggle = 1;
            for (int j = 0; j < n - 1; ++j) // каждый столбец
            {
                double colMax = Math.Abs(result[j][j]); // Наибольшее значение в столбце j
                int pRow = j;
                for (int i = j + 1; i < n; ++i)
                {
                    if (result[i][j] > colMax)
                    {
                        colMax = result[i][j];
                        pRow = i;
                    }
                }
                if (pRow != j) // перестановка строк
                {
                    double[] rowPtr = result[pRow];
                    result[pRow] = result[j];
                    result[j] = rowPtr;
                    int tmp = perm[pRow]; // Меняем информацию о перестановке
                    perm[pRow] = perm[j];
                    perm[j] = tmp;
                    toggle = -toggle; // переключатель перестановки строк
                }
                if (Math.Abs(result[j][j]) < 1.0E-20)
                    return null;
                for (int i = j + 1; i < n; ++i)
                {
                    result[i][j] /= result[j][j];
                    for (int k = j + 1; k < n; ++k)
                        result[i][k] -= result[i][j] * result[j][k];
                }
            } // основной цикл по столбцу j
            double[,] resul = new double[n, n];
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    resul[i, j] = result[i][j];

            return resul;
        }
        public static double[] HelperSolve(double[,] luMatrix, double[] b)
        {
            // Решаем luMatrix * x = b
            int n = 6;
            double[] x = new double[n];
            b.CopyTo(x, 0);
            for (int i = 1; i < n; ++i)
            {
                double sum = x[i];
                for (int j = 0; j < i; ++j)
                    sum -= luMatrix[i, j] * x[j];
                x[i] = sum;
            }
            x[n - 1] /= luMatrix[n - 1, n - 1];
            for (int i = n - 2; i >= 0; --i)
            {
                double sum = x[i];
                for (int j = i + 1; j < n; ++j)
                    sum -= luMatrix[i, j] * x[j];
                x[i] = sum / luMatrix[i, i];
            }
            return x;
        }
        public static double[,] MatrixInverse(double[,] matrix)
        {
            int n = 6;
            double[,] result = matrix;
            int[] perm;
            int toggle;
            double[,] lum = MatrixDecompose(matrix, out perm, out toggle);
            if (lum == null)
                throw new Exception("Невозможно найти обратную матрицу");
            double[] b = new double[n];
            for (int i = 0; i < n; ++i)
            {
                for (int j = 0; j < n; ++j)
                {
                    if (i == perm[j])
                        b[j] = 1.0;
                    else
                        b[j] = 0.0;
                }
                double[] x = HelperSolve(lum, b);
                for (int j = 0; j < n; ++j)
                    result[j, i] = x[j];
            }
            return result;
        }

        //Aabcd*Bcd
        public static double[,] MulSvert(double[,,,] tenzor1, double[,] tenzor2)
        {
            double[,] mulsvertenzor = new double[3, 3];
            Console.WriteLine("\n\nУмножение со свёрткой:\n");
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    for (int k = 0; k < 3; k++)
                        for (int l = 0; l < 3; l++)
                            mulsvertenzor[i, j] += tenzor1[i, j, k, l] * tenzor2[k, l];
                    Console.Write(mulsvertenzor[i, j] + " ");
                }
                Console.WriteLine();
            }
            return mulsvertenzor;
        }
    }

    class Program
    {
        static void Main(string[] args)
        {
            //Пункт 0 (ввод данных)
            Console.WriteLine("Введите значение скорости продольных волн Vp: ");
            double Vp = double.Parse(Console.ReadLine());
            Console.WriteLine("\nВведите значение скорости поперечных волн Vs: ");
            double Vs = double.Parse(Console.ReadLine());
            Console.WriteLine("\nВведите значение плотности образца Ro: ");
            double Ro = double.Parse(Console.ReadLine());
            Console.WriteLine("\nВведите значения концентраций минерала f1: ");
            int f1 = int.Parse(Console.ReadLine());
            Console.WriteLine("\nи жидкости f2: ");
            int f2 = int.Parse(Console.ReadLine());
            Console.WriteLine("\nВведите значение модуля сдвига минерала мю1: ");
            int Mum1 = int.Parse(Console.ReadLine());
            Console.WriteLine("\nВведите значение среднего давления жидкости в порах минерала P: ");
            int Pwater = int.Parse(Console.ReadLine());
            const int Mum2 = 0;
            //______________________________________________________________________________________________________________________________________________________________________________________________________________________

            //Пункт 1 (расчет Sigma(0g))
            double[,] sigma0g = new double[3, 3];
            sigma0g = Sigma.Calculation();
            Console.WriteLine("\n\nsigma0g:\n");
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                    Console.Write(sigma0g[i, j] + " ");
                Console.WriteLine();
            }

            //______________________________________________________________________________________________________________________________________________________________________________________________________________________

            //Пункт 2 (нахождение С(m))
            double Mum = Math.Pow(Vs, 2) * Ro;
            double Lyamm = Math.Pow(Vp, 2) * Ro - 2 * Mum;

            double[,] Cm = new double[6, 6];
            Console.WriteLine("\n\nCm: \n");
            for (int i = 0; i < 6; i++)
            {
                if (i < 3)
                    Cm[i, i] = Lyamm + 2 * Mum;
                else
                    Cm[i, i] = Mum;

                for (int j = 0; j < 6; j++)
                {
                    if ((i != j) && (i < 3) && (j < 3))
                        Cm[i, j] = Lyamm;
                    else if (i != j)
                        Cm[i, j] = 0;
                    Console.Write(Cm[i, j] + "\t");
                }
                Console.WriteLine();
            }

            //______________________________________________________________________________________________________________________________________________________________________________________________________________________

            //Пункт 3 (преобразование в Cm ijkl)
            double[,,,] CmTenzor = new double[3, 3, 3, 3];
            //1-11      2-22       3-33       4-23 32        5-31 13        6-12 21
            for (int i = 0; i < 6; i++)
                for (int j = 0; j < 6; j++)
                {
                    if ((i < 3) && (j < 3))                                             //11, 12, 13, 21, 22, 23, 31, 32, 33
                        CmTenzor[i, i, j, j] = Cm[i, j];
                    else if (i < 3)                                                     //14, 15, 16, 24, 25, 26, 34, 35, 36
                    {
                        if (j == 3)                                                     //14, 24, 34
                            CmTenzor[i, i, 1, 2] = CmTenzor[i, i, 2, 1] = Cm[i, j];
                        else if (j == 4)                                                //15, 25, 35
                            CmTenzor[i, i, 2, 0] = CmTenzor[i, i, 0, 2] = Cm[i, j];
                        else                                                            //16, 26, 36
                            CmTenzor[i, i, 0, 1] = CmTenzor[i, i, 1, 0] = Cm[i, j];
                    }
                    else if (j < 3)                                                     //41, 42, 43, 51, 52, 53, 61, 62, 63
                    {
                        if (i == 3)                                                     //41, 42, 43
                            CmTenzor[1, 2, j, j] = CmTenzor[2, 1, j, j] = Cm[i, j];
                        else if (i == 4)                                                //51, 52, 53
                            CmTenzor[2, 0, j, j] = CmTenzor[0, 2, j, j] = Cm[i, j];
                        else                                                            //61, 62, 63
                            CmTenzor[0, 1, j, j] = CmTenzor[1, 0, j, j] = Cm[i, j];
                    }
                    else if (i == 3)                                                    //44, 45, 46
                    {
                        if (j == 3)                                                     //44
                            CmTenzor[1, 2, 1, 2] = CmTenzor[2, 1, 2, 1] = CmTenzor[1, 2, 2, 1] = CmTenzor[2, 1, 1, 2] = Cm[i, j];
                        else if (j == 4)                                                //45
                            CmTenzor[1, 2, 0, 2] = CmTenzor[2, 1, 2, 0] = CmTenzor[1, 2, 2, 0] = CmTenzor[2, 1, 0, 2] = Cm[i, j];
                        else                                                            //46
                            CmTenzor[1, 2, 0, 1] = CmTenzor[2, 1, 1, 0] = CmTenzor[1, 2, 1, 0] = CmTenzor[2, 1, 0, 1] = Cm[i, j];
                    }
                    else if (i == 4)                                                    //54, 55, 56
                    {
                        if (j == 3)                                                     //54
                            CmTenzor[0, 2, 1, 2] = CmTenzor[2, 0, 2, 1] = CmTenzor[0, 2, 2, 1] = CmTenzor[2, 0, 1, 2] = Cm[i, j];
                        else if (j == 4)                                                //55
                            CmTenzor[0, 2, 0, 2] = CmTenzor[2, 0, 2, 0] = CmTenzor[0, 2, 2, 0] = CmTenzor[2, 0, 0, 2] = Cm[i, j];
                        else                                                            //56
                            CmTenzor[0, 2, 0, 1] = CmTenzor[2, 0, 1, 0] = CmTenzor[0, 2, 1, 0] = CmTenzor[2, 0, 0, 1] = Cm[i, j];
                    }
                    else                                                                //64, 65, 66
                    {
                        if (j == 3)                                                     //64
                            CmTenzor[1, 0, 1, 2] = CmTenzor[0, 1, 2, 1] = CmTenzor[1, 0, 2, 1] = CmTenzor[0, 1, 1, 2] = Cm[i, j];
                        else if (j == 4)                                                //65
                            CmTenzor[1, 0, 0, 2] = CmTenzor[0, 1, 2, 0] = CmTenzor[1, 0, 2, 0] = CmTenzor[0, 1, 0, 2] = Cm[i, j];
                        else                                                            //66
                            CmTenzor[1, 0, 0, 1] = CmTenzor[0, 1, 1, 0] = CmTenzor[1, 0, 1, 0] = CmTenzor[0, 1, 0, 1] = Cm[i, j];
                    }
                }

            Console.WriteLine("\n\nCijkl:\n");
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                    for (int k = 0; k < 3; k++)
                        for (int l = 0; l < 3; l++)
                            Console.Write(CmTenzor[i, j, k, l] + "  ");
                Console.WriteLine();
            }
            Console.WriteLine();
            //______________________________________________________________________________________________________________________________________________________________________________________________________________________

            //Пункт 4 (нахождение С(mwp))

            // Voight
            int Muv = f1 * Mum1 + f2 * Mum2;
            int Km1 = (int)Lyamm + 2 / 3 * (int)Mum;
            int Km2 = Pwater;
            int Kmv = f1 * Km1 + f2 * Km2;

            // Roice
            int Mur = Mum2;
            int Kmr = (int)Math.Pow((f1 / Km1 + f2 / Km2), (-1));

            // границы 
            var rnd = new Random();

            // Kmr <= Kmwp <= Kmv
            double Kmwp = (double)((rnd.Next(Kmr, Kmv) / 1.0f));

            // Mur <= Mumwp <= Muv
            double Mumwp = (double)((rnd.Next(Mur, Muv) / 1.0f));
            Console.Write("\n\nK = {0}, мю = {1}\n", Kmwp, Mumwp);

            double Lyammwp = Kmwp - 2 / 3 * Mumwp;

            double[,] Cmwp = new double[6, 6];
            Console.WriteLine("\n\nCmwp:\n");
            for (int i = 0; i < 6; i++)
            {
                if (i < 3)
                    Cmwp[i, i] = Lyammwp + 2 * Muv;
                else
                    Cmwp[i, i] = Muv;

                for (int j = 0; j < 6; j++)
                {
                    if ((i != j) && (i < 3) && (j < 3))
                        Cmwp[i, j] = Lyammwp;
                    else if (i != j)
                        Cmwp[i, j] = 0;
                    Console.Write(Cmwp[i, j] + "\t");
                }
                Console.WriteLine();
            }

            //______________________________________________________________________________________________________________________________________________________________________________________________________________________

            //Пункт 5 (нахождение S(mwp))
            double[,] Smwp = new double[6, 6];
            Smwp = Matrix.MatrixInverse(Cmwp);
            Console.WriteLine("\n\nSmwp:\n");
            for (int i = 0; i < 6; i++)
            {
                for (int j = 0; j < 6; j++)
                    Console.Write("{0,6:f2}   ", Smwp[i, j]);
                Console.WriteLine();
            }

            //______________________________________________________________________________________________________________________________________________________________________________________________________________________

            //Пункт 6 (нахождение S(mwp)ijkl)
            double[,,,] SmwpTenzor = new double[3, 3, 3, 3];
            //1-11      2-22       3-33       4-23 32        5-31 13        6-12 21
            for (int i = 0; i < 6; i++)
                for (int j = 0; j < 6; j++)
                {
                    if ((i < 3) && (j < 3))                                             //11, 12, 13, 21, 22, 23, 31, 32, 33
                        SmwpTenzor[i, i, j, j] = Smwp[i, j];
                    else if (i < 3)                                                     //14, 15, 16, 24, 25, 26, 34, 35, 36
                    {
                        if (j == 3)                                                     //14, 24, 34
                            SmwpTenzor[i, i, 1, 2] = SmwpTenzor[i, i, 2, 1] = Smwp[i, j] / 2;
                        else if (j == 4)                                                //15, 25, 35
                            SmwpTenzor[i, i, 2, 0] = SmwpTenzor[i, i, 0, 2] = Smwp[i, j] / 2;
                        else                                                            //16, 26, 36
                            SmwpTenzor[i, i, 0, 1] = SmwpTenzor[i, i, 1, 0] = Smwp[i, j] / 2;
                    }
                    else if (j < 3)                                                     //41, 42, 43, 51, 52, 53, 61, 62, 63
                    {
                        if (i == 3)                                                     //41, 42, 43
                            SmwpTenzor[1, 2, j, j] = SmwpTenzor[2, 1, j, j] = Smwp[i, j] / 2;
                        else if (i == 4)                                                //51, 52, 53
                            SmwpTenzor[2, 0, j, j] = SmwpTenzor[0, 2, j, j] = Smwp[i, j] / 2;
                        else                                                            //61, 62, 63
                            SmwpTenzor[0, 1, j, j] = SmwpTenzor[1, 0, j, j] = Smwp[i, j] / 2;
                    }
                    else if (i == 3)                                                    //44, 45, 46
                    {
                        if (j == 3)                                                     //44
                            SmwpTenzor[1, 2, 1, 2] = SmwpTenzor[2, 1, 2, 1] = SmwpTenzor[1, 2, 2, 1] = SmwpTenzor[2, 1, 1, 2] = Smwp[i, j] / 4;
                        else if (j == 4)                                                //45
                            SmwpTenzor[1, 2, 0, 2] = SmwpTenzor[2, 1, 2, 0] = SmwpTenzor[1, 2, 2, 0] = SmwpTenzor[2, 1, 0, 2] = Smwp[i, j] / 4;
                        else                                                            //46
                            SmwpTenzor[1, 2, 0, 1] = SmwpTenzor[2, 1, 1, 0] = SmwpTenzor[1, 2, 1, 0] = SmwpTenzor[2, 1, 0, 1] = Smwp[i, j] / 4;
                    }
                    else if (i == 4)                                                    //54, 55, 56
                    {
                        if (j == 3)                                                     //54
                            SmwpTenzor[0, 2, 1, 2] = SmwpTenzor[2, 0, 2, 1] = SmwpTenzor[0, 2, 2, 1] = SmwpTenzor[2, 0, 1, 2] = Smwp[i, j] / 4;
                        else if (j == 4)                                                //55
                            SmwpTenzor[0, 2, 0, 2] = SmwpTenzor[2, 0, 2, 0] = SmwpTenzor[0, 2, 2, 0] = SmwpTenzor[2, 0, 0, 2] = Smwp[i, j] / 4;
                        else                                                            //56
                            SmwpTenzor[0, 2, 0, 1] = SmwpTenzor[2, 0, 1, 0] = SmwpTenzor[0, 2, 1, 0] = SmwpTenzor[2, 0, 0, 1] = Smwp[i, j] / 4;
                    }
                    else                                                                //64, 65, 66
                    {
                        if (j == 3)                                                     //64
                            SmwpTenzor[1, 0, 1, 2] = SmwpTenzor[0, 1, 2, 1] = SmwpTenzor[1, 0, 2, 1] = SmwpTenzor[0, 1, 1, 2] = Smwp[i, j] / 4;
                        else if (j == 4)                                                //65
                            SmwpTenzor[1, 0, 0, 2] = SmwpTenzor[0, 1, 2, 0] = SmwpTenzor[1, 0, 2, 0] = SmwpTenzor[0, 1, 0, 2] = Smwp[i, j] / 4;
                        else                                                            //66
                            SmwpTenzor[1, 0, 0, 1] = SmwpTenzor[0, 1, 1, 0] = SmwpTenzor[1, 0, 1, 0] = SmwpTenzor[0, 1, 0, 1] = Smwp[i, j] / 4;
                    }
                }

            Console.WriteLine("\n\nSijkl:\n");
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                    for (int k = 0; k < 3; k++)
                        for (int l = 0; l < 3; l++)
                            Console.Write("{0,6:f2}  ", SmwpTenzor[i, j, k, l]);
                Console.WriteLine();
            }
            Console.WriteLine();


            //______________________________________________________________________________________________________________________________________________________________________________________________________________________

            //Пункт 7 (нахождение Sigma(m)ij)
            double[,] SigmamTenzor = new double[3, 3];
            // Sigmam mn = Cm mnij * (Smwp ijkl * Sigma0g kl)

            double[,] prom = Matrix.MulSvert(SmwpTenzor, sigma0g);
            SigmamTenzor = Matrix.MulSvert(CmTenzor, prom);

            Console.WriteLine("\n\nSigma(m):\n");
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                    for (int k = 0; k < 3; k++)
                        for (int l = 0; l < 3; l++)
                            Console.Write("{0,6:f2}  ", SigmamTenzor[i, j]);
                Console.WriteLine();
            }
            Console.WriteLine();

            //______________________________________________________________________________________________________________________________________________________________________________________________________________________

            Console.ReadKey();
        }

    }

    class Sugl
    {
        public double plotnost = 1750;
        public double koefPD = 0.35;
    }
    class Supes
    {
        public double plotnost = 1750;
        public double koefPD = 0.3;
    }
    class Lyos
    {
        public double plotnost = 2200;
        public double koefPD = 0.25;
    }
    public class Sigma
    {
        public static double[,] Calculation()
        {
            const double g = 9.81;
            Sugl suglinki = new Sugl();
            Supes supesi = new Supes();
            Lyos lyosi = new Lyos();

            //Console.WriteLine(" _______________ _______________ _______________ _________________________________ ______________ ");
            //Console.WriteLine("| Отдел \t| Ярус \t\t| Индекс \t| Характеристики пород \t\t  | Мощность в м |");
            //Console.WriteLine("|_______________|_______________|_______________|_________________________________|______________|");
            //Console.WriteLine("| Современный \t| — \t \t| (Q4)^2 (1) \t| Суглинки, супеси \t\t  | 5-6 \t |");
            //Console.WriteLine("| Современный \t| — \t \t| (Q4)^1 (2) \t| Суглинки, супеси \t\t  | до 10 \t |");
            //Console.WriteLine("| Верхний \t| Хвалынский \t| Q3hv^2 (3) \t| Лёссы, супеси \t\t  | 20+-1 \t |");
            //Console.WriteLine("|_______________|_______________|_______________|_________________________________|______________|\n");

            //Console.WriteLine("\nВведите мощность для всех указанных слоёв:");
            double[] power = new double[3];
            power[0] = 5;
            power[1] = 9;
            power[2] = 20;
            //for (int i = 0; i < 3; i++)
            //    power[i] = double.Parse(Console.ReadLine());

            //Console.WriteLine("\nВведите наименование горной породы строчными буквами:");
            string nameGP = "супеси"; /*= Console.ReadLine();*/
            //Console.WriteLine("\nЖелаете ли вы воспользоваться имеющимися данными?");
            //Console.WriteLine("Введите 1, если \"да\", и 0, если \"нет\"");
            int ch = 1; /*= int.Parse(Console.ReadLine());*/

            if (ch == 0)
            {
                Console.WriteLine("\nВведите плотность для породы суглинки:");
                suglinki.plotnost = double.Parse(Console.ReadLine());

                Console.WriteLine("\nВведите плотность для породы лёссы:");
                lyosi.plotnost = double.Parse(Console.ReadLine());

                Console.WriteLine("\nВведите плотность для породы супеси:");
                supesi.plotnost = double.Parse(Console.ReadLine());

                if (nameGP == "суглинки")
                {
                    Console.WriteLine("\nВведите коэффициент Пуассона для породы суглинки:");
                    suglinki.koefPD = double.Parse(Console.ReadLine());
                }
                else if (nameGP == "лёссы")
                {
                    Console.WriteLine("\nВведите коэффициент Пуассона для породы лёссы:");
                    lyosi.koefPD = double.Parse(Console.ReadLine());
                }
                else if (nameGP == "супеси")
                {
                    Console.WriteLine("\nВведите коэффициент Пуассона для породы супеси:");
                    supesi.koefPD = double.Parse(Console.ReadLine());
                }
                else
                    Console.WriteLine("Данной породы нет в пласте");
            }
            else
            {
                Console.WriteLine("Плотность породы суглинки:  " + suglinki.plotnost);
                Console.WriteLine("\nПлотность породы супеси:  " + supesi.plotnost);
                Console.WriteLine("\nПлотность породы лёссы:  " + lyosi.plotnost);

                if (nameGP == "суглинки")
                    Console.WriteLine("\nКоэффициент Пуассона для породы суглинки:  " + suglinki.koefPD);
                else if (nameGP == "супеси")
                    Console.WriteLine("\nКоэффициент Пуассона для породы супеси:  " + supesi.koefPD);
                else if (nameGP == "лёссы")
                    Console.WriteLine("\nКоэффициент Пуассона для породы лёссы:  " + lyosi.koefPD);
                else
                    Console.WriteLine("\nДанной породы нет в пласте");
            }

            double[] gamma = new double[3];
            gamma[0] = (suglinki.plotnost + supesi.plotnost) / 2 * g;
            gamma[1] = (suglinki.plotnost + supesi.plotnost) / 2 * g;
            gamma[2] = (lyosi.plotnost + supesi.plotnost) / 2 * g;

            double sigma3 = 0;
            for (int i = 0; i < 3; i++)
                sigma3 += gamma[i] * power[i];

            double koefBR = 0;
            if (nameGP == "лёссы")
            {
                koefBR = lyosi.koefPD / (1 - lyosi.koefPD);
                Console.WriteLine("\nКоэффициент бокового распора:  " + koefBR);
            }
            else if (nameGP == "супеси")
            {
                koefBR = supesi.koefPD / (1 - supesi.koefPD);
                Console.WriteLine("\nКоэффициент бокового распора:  " + koefBR);
            }
            else
                Console.WriteLine("\nДанной породы нет в таблице");

            double sigma1, sigma2;
            sigma1 = koefBR * sigma3;
            sigma2 = sigma1;
            double[,] sigma = new double[3, 3] { { sigma1, 0, 0 }, { 0, sigma2, 0 }, { 0, 0, sigma3 } };

            //Console.WriteLine("\nТектоническая составляющая равна: \n");
            //for (int i = 0; i < 3; i++)
            //{
            //    for (int j = 0; j < 3; j++)
            //        Console.Write(sigma[i, j] + "\t");
            //    Console.WriteLine();
            //}

            return sigma;
        }
    }
}