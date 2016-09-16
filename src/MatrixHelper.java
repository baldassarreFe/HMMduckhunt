import java.io.IOException;
import java.io.InputStreamReader;
import java.text.NumberFormat;
import java.util.Arrays;
import java.util.Random;
import java.util.StringTokenizer;

public class MatrixHelper {

	public static double[][] transpose(double[][] matrix) {
		int rows = matrix.length;
		int columns = matrix[0].length;
		double[][] result = new double[columns][rows];
		for (int i = 0; i < rows; i++)
			for (int j = 0; j < columns; j++)
				result[j][i] = matrix[i][j];
		return result;
	}

	public double[][] sum(double[][] first, double[][] second) {
		int rows = first.length;
		int columns = first[0].length;
		double[][] result = new double[rows][columns];
		for (int row = 0; row < result.length; row++) {
			for (int column = 0; column < result[0].length; column++) {
				result[row][column] = first[row][column] + second[row][column];
			}
		}
		return result;
	}

	public double[][] elementMultiply(double[][] first, double[][] second) {
		int rows = first.length;
		int columns = first[0].length;
		double[][] result = new double[rows][columns];
		for (int row = 0; row < result.length; row++) {
			for (int column = 0; column < result[0].length; column++) {
				result[row][column] = first[row][column] * second[row][column];
			}
		}
		return result;
	}

	public static double[][] multiply(double[][] first, double[][] second) {
		int firstRows = first.length;
		int firstColumns = first[0].length;
		int secondRows = second.length;
		int secondColumns = second[0].length;

		check(firstColumns == secondRows);

		double[][] result = new double[firstRows][secondColumns];
		for (int firstRow = 0; firstRow < firstRows; firstRow++) {
			for (int secondColumn = 0; secondColumn < secondColumns; secondColumn++) {
				result[firstRow][secondColumn] = multiply(getRow(first, firstRow), getColumn(second, secondColumn));
			}
		}

		return result;
	}

	public static double multiply(double[] first, double[] second) {
		check(first.length == second.length);
		double t = 0;
		for (int i = 0; i < first.length; i++) {
			t += first[i] * second[i];
		}
		return t;
	}

	private static void check(boolean condition) {
		if (!condition)
			throw new IllegalArgumentException();
	}

	public static double[][] power(double[][] matrix, int exponent) {
		check(matrix.length == matrix[0].length);
		double[][] result = identity(matrix.length);
		while (exponent-- > 0) {
			result = multiply(result, matrix);
		}
		return result;
	}

	public static double[] getRow(double[][] matrix, int row) {
		double[] result = Arrays.copyOf(matrix[row], matrix[0].length);
		return result;
	}

	public static double[] getColumn(double[][] matrix, int column) {
		double[] result = new double[matrix.length];
		for (int i = 0; i < matrix.length; i++) {
			result[i] = matrix[i][column];
		}
		return result;
	}

	public static int indexOfHigherElement(double[] array) {
		int index = 0;
		double higher = array[0];
		for (int i = 0; i < array.length; i++) {
			if (array[i] > higher) {
				index = i;
			}
		}
		return index;
	}

	public static String toString(double[][] matrix) {
		StringBuilder sb = new StringBuilder();
		for (int row = 0; row < matrix.length; row++) {
			sb.append("[");
			for (int column = 0; column < matrix[0].length; column++) {
				sb.append("\t" + matrix[row][column]);
			}
			sb.append("\t]\n");
		}
		return sb.toString();
	}

	public static double[][] identity(int size) {
		double[][] result = new double[size][size];
		for (int i = 0; i < size; i++) {
			result[i][i] = 1;
		}
		return result;
	}
	
	public static double[][] newAlmostDiagonalMatrix(int size) {
		double[][] matrix = new double[size][size];
		almostDiagonalMatrix(matrix);
		return matrix;
	}

	public static void almostDiagonalMatrix(double[][] matrix) {
		Random r = new Random();
		int rows = matrix.length;
		int columns = matrix[0].length;
		for (int row = 0; row < rows; row++) {
			double factor = 0;
			for (int column = 0; column < columns; column++) {
				matrix[row][column] = 1.0 / columns * ((row == column ? 8 : 1) + r.nextDouble());
				factor += matrix[row][column];
			}
			for (int column = 0; column < columns; column++) {
				matrix[row][column] /= factor;
			}
		}
	}

	public static double[][] newRowStochasticMatrix(int rows, int columns) {
		double[][] matrix = new double[rows][columns];
		rowStochasticMatrix(matrix);
		return matrix;
	}

	public static double[] newRowStochasticArray(int columns) {
		double[] array = new double[columns];
		rowStochasticArray(array);
		return array;
	}

	public static void rowStochasticMatrix(double[][] matrix) {
		int rows = matrix.length;
		for (int row = 0; row < rows; row++) {
			rowStochasticArray(matrix[row]);
		}
	}

	public static void rowStochasticArray(double[] array) {
		Random r = new Random();
		double factor = 0;
		for (int column = 0; column < array.length; column++) {
			array[column] = r.nextDouble();
			factor += array[column];
		}
		for (int column = 0; column < array.length; column++) {
			array[column] /= factor;
		}
	}

	public static int argMax(double[] array) {
		int bestIndex = 0;
		for (int i = 0; i < array.length; i++) {
			if (array[i] > array[bestIndex])
				bestIndex = i;
		}
		return bestIndex;
	}

	public static double[][] parseMatrix(String line) {
		StringTokenizer st = new StringTokenizer(line, " ");
		int rows = Integer.parseInt(st.nextToken());
		int columns = Integer.parseInt(st.nextToken());
		double[][] result = new double[rows][columns];
		for (int row = 0; row < rows; row++) {
			for (int column = 0; column < columns; column++) {
				result[row][column] = Double.parseDouble(st.nextToken());
			}
		}
		return result;
	}

	public static int[] parseIntArray(String line) {
		StringTokenizer st = new StringTokenizer(line, " ");
		int columns = Integer.parseInt(st.nextToken());
		int[] result = new int[columns];
		for (int column = 0; column < columns; column++) {
			result[column] = Integer.parseInt(st.nextToken());
		}
		return result;
	}

	public static int[] parseIntArray(InputStreamReader input) throws IOException {
		int lenght = readInt(input);
		int[] result = new int[lenght];
		for (int i = 0; i < result.length; i++) {
			result[i] = readInt(input);
		}
		return result;
	}

	private static int readInt(InputStreamReader input) throws IOException {
		int c;
		StringBuilder sb = new StringBuilder();
		do {
			c = input.read();
			sb.append((char) c);
		} while (c >= '0' && c <= '9');
		return Integer.parseInt(sb.toString().trim());
	}

	public static String stringifyMatrix(double[][] matrix) {
		StringBuilder sb = new StringBuilder();
		sb.append(matrix.length);
		sb.append(" ");
		sb.append(matrix[0].length);
		for (int row = 0; row < matrix.length; row++) {
			for (int column = 0; column < matrix[0].length; column++) {
				sb.append(" ");
				sb.append(matrix[row][column]);
			}
		}
		return sb.toString();
	}

	public static String stringifyArray(double[] array, boolean printSize) {
		StringBuilder sb = new StringBuilder();
		if (printSize) {
			sb.append("1 ");
			sb.append(array.length);
		}
		for (int column = 0; column < array.length; column++) {
			sb.append(" ");
			sb.append(array[column]);
		}
		return sb.toString();
	}

	public static String stringifyArray(int[] array) {
		return stringifyArray(array, true);
	}

	public static String stringifyArray(int[] array, boolean printSize) {
		StringBuilder sb = new StringBuilder();
		if (printSize) {
			sb.append("1 ");
			sb.append(array.length);
		}
		for (int column = 0; column < array.length; column++) {
			sb.append(" ");
			sb.append(array[column]);
		}
		return sb.toString();
	}

	public static String stringifyArray(double[] array) {
		return stringifyArray(array, true);
	}

	public static String toPrettyString(double[][] matrix) {
		StringBuilder sb = new StringBuilder();
		for (int row = 0; row < matrix.length; row++) {
			sb.append(toPrettyString(matrix[row]));
		}
		return sb.toString();
	}

	public static String toPrettyString(double[] array) {
		NumberFormat nf = NumberFormat.getNumberInstance();
		nf.setMinimumFractionDigits(4);
		nf.setMaximumFractionDigits(4);
		StringBuilder sb = new StringBuilder();
		sb.append("[");
		for (int column = 0; column < array.length; column++) {
			String number;
			if (Double.isNaN(array[column]))
				number = "NaN";
			else
				number = nf.format(array[column]);

			sb.append("\t" + number);
		}
		sb.append("\t]\n");
		return sb.toString();
	}

	public static double distance(double[][] first, double[][] second) {
		check(first.length == second.length && first[0].length == second[0].length);
		double dist = 0;
		for (int i = 0; i < first.length; i++) {
			for (int j = 0; j < first[0].length; j++) {
				dist += Math.abs(first[i][j] - second[i][j]);
			}
		}
		return dist;
	}

	public static double distance(double[] first, double[] second) {
		check(first.length == second.length);
		double dist = 0;
		for (int i = 0; i < first.length; i++) {
			dist += Math.pow((first[i] - second[i]), 2);
		}
		return Math.sqrt(dist);
	}

	public static void copy(double[][] from, double[][] to) {
		for (int row = 0; row < from.length; row++) {
			for (int column = 0; column < from[0].length; column++) {
				to[row][column] = from[row][column];
			}
		}
	}

	public static void copy(double[] from, double[] to) {
		for (int column = 0; column < from.length; column++) {
			to[column] = from[column];
		}
	}

	public static double[] normalize(double[] array) {
		double[] result = new double[array.length];
		double factor = 0;
		for (int column = 0; column < array.length; column++) {
			factor += array[column];
		}
		if (factor != 0)
			for (int column = 0; column < array.length; column++) {
				result[column] = array[column]/factor;
			}
		return result;
	}
}
