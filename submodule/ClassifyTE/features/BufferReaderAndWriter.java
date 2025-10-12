
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;

public class BufferReaderAndWriter {
	
	// this method creates the buffer reader and return it back to the calling method
	public static BufferedReader getReader(File file) {

		FileInputStream fstream;
		try {
			FileInputStream fis = new FileInputStream(file);
			InputStreamReader isr = new InputStreamReader(fis);
			BufferedReader br = new BufferedReader(isr);
			return br;

		} catch (Exception e) {
			return null;
		}
	}

	// this method creates the buffer writer and returns it back to the calling method
	
	public static BufferedWriter getWriter(File file) {

		FileOutputStream fstream;
		try {
			FileOutputStream fos = new FileOutputStream(file);
			OutputStreamWriter osw = new OutputStreamWriter(fos, "UTF-8");
			BufferedWriter bufferedWriter = new BufferedWriter(osw);
			return bufferedWriter;

		} catch (Exception e) {
			return null;
		}
	}

}
