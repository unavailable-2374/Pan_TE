import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;


public class KmersFeaturesCollector {

	public static void main(String[] args) throws IOException {
		
		File currentDir = new File(new File(".").getAbsolutePath());
		String curr_dir_path = currentDir.getCanonicalPath();
		
		String listfilePath = curr_dir_path+"/list.txt";
		//System.out.println(listfilePath);
		BufferedReader list_rd = BufferReaderAndWriter.getReader(new File(listfilePath));
		String id = "";
		String output_feature_file_path = curr_dir_path+"/feature_file.csv";
		FileWriter output_feature_file = new FileWriter(output_feature_file_path, true);
		BufferedWriter feature_bw = new BufferedWriter(output_feature_file);
		PrintWriter feature_wr = new PrintWriter(feature_bw);
		
		// initialize the hasmap and insert keys
		Map<String, Integer> feature_map = new HashMap<String, Integer>();	
		List<String> neuc_type = new ArrayList<String>();
		neuc_type.add("A");
		neuc_type.add("C");
		neuc_type.add("T");
		neuc_type.add("G");
		
		StringBuilder header = new StringBuilder();
		



		// generate two mer types
		for(int i = 0; i < 4; i++){
			for(int j = 0; j < 4; j++){
				
				String two_mer_type = neuc_type.get(i)+neuc_type.get(j);
				feature_map.put(two_mer_type, 0);
				header.append(two_mer_type+",");
				
			}			
			
		}
		
		// generate three mer types
		for(int i = 0; i < 4; i++){
			for(int j = 0; j < 4; j++){
				for(int k = 0; k < 4; k++){
					
					String three_mer_type = neuc_type.get(i)+neuc_type.get(j)+neuc_type.get(k);
					feature_map.put(three_mer_type, 0);
					header.append(three_mer_type+",");
					
				}
				
				
			}			
			
		}
		
		// generate four mer types
		for(int i = 0; i < 4; i++){
			for(int j = 0; j < 4; j++){
				for(int k = 0; k < 4; k++){
					for(int l = 0; l < 4; l++){
						
						String four_mer_type = neuc_type.get(i)+neuc_type.get(j)+neuc_type.get(k)+neuc_type.get(l);
						feature_map.put(four_mer_type, 0);
						header.append(four_mer_type+",");
														
					}
					
				}
				
				
			}			
			
		}
		

		String header_types = header.toString();
		String header_final = header_types.substring(0, header_types.length()-1);
		feature_wr.println(header_final);
		feature_wr.flush();
		
		while((id = list_rd.readLine())!=null){
			
			String two_mer_file = curr_dir_path+"/kanalyze-2.0.0/output_data/2mer/"+id+".txt";
//			System.out.println(two_mer_file);
			String three_mer_file = curr_dir_path+"/kanalyze-2.0.0/output_data/3mer/"+id+".txt";
			String four_mer_file = curr_dir_path+"/kanalyze-2.0.0/output_data/4mer/"+id+".txt";
			
			
			for (Map.Entry<String, Integer> entry : feature_map.entrySet()) {
		    
				feature_map.put(entry.getKey(), 0);
		    
			}
			
			parseInputAndAppendToOutput(two_mer_file, three_mer_file, four_mer_file, feature_wr, feature_map, neuc_type);
			
			
		}
		
		list_rd.close();
		feature_wr.close();
		
		
	}
	
	
	public static void parseInputAndAppendToOutput(String two_mer_file, String three_mer_file, String four_mer_file, PrintWriter feature_wr, Map<String, Integer> feature_map, List<String> neuc_type) throws IOException{
		
		
		BufferedReader two_mer_rd = BufferReaderAndWriter.getReader(new File(two_mer_file));
		String two_mer_line;
		BufferedReader three_mer_rd = BufferReaderAndWriter.getReader(new File(three_mer_file));
		String three_mer_line;
		BufferedReader four_mer_rd = BufferReaderAndWriter.getReader(new File(four_mer_file));
		String four_mer_line;
		

		while((two_mer_line = two_mer_rd.readLine())!=null){
			
			String[] tokens = two_mer_line.split("\t");
			String key = tokens[0];
			Integer value = Integer.parseInt(tokens[1]);
			
			feature_map.put(key, value);
			
		}
		
		while((three_mer_line = three_mer_rd.readLine())!=null){
			
			String[] tokens = three_mer_line.split("\t");
			String key = tokens[0];
			Integer value = Integer.parseInt(tokens[1]);
			
			feature_map.put(key, value);
			
		}
		
		while((four_mer_line = four_mer_rd.readLine())!=null){
			
			String[] tokens = four_mer_line.split("\t");
			String key = tokens[0];
			Integer value = Integer.parseInt(tokens[1]);
			
			feature_map.put(key, value);
			
		}
		
		
		
		two_mer_rd.close();
		three_mer_rd.close();
		four_mer_rd.close();
		
		
		StringBuilder featues = new StringBuilder();
		
		


		// get two mer values
		for(int i = 0; i < 4; i++){
			for(int j = 0; j < 4; j++){
				
				String two_mer_type = neuc_type.get(i)+neuc_type.get(j);
				featues.append(feature_map.get(two_mer_type)+",");
				
			}			
			
		}
		
		// generate three mer types
		for(int i = 0; i < 4; i++){
			for(int j = 0; j < 4; j++){
				for(int k = 0; k < 4; k++){
					
					String three_mer_type = neuc_type.get(i)+neuc_type.get(j)+neuc_type.get(k);
					featues.append(feature_map.get(three_mer_type)+",");
					
				}
				
				
			}			
			
		}
		
		// generate four mer types
		for(int i = 0; i < 4; i++){
			for(int j = 0; j < 4; j++){
				for(int k = 0; k < 4; k++){
					for(int l = 0; l < 4; l++){
						
						String four_mer_type = neuc_type.get(i)+neuc_type.get(j)+neuc_type.get(k)+neuc_type.get(l);
						featues.append(feature_map.get(four_mer_type)+",");
														
					}
					
				}
				
				
			}			
			
		}
		
		
		
		
		String feature_string = featues.toString();
		String feature_final = feature_string.substring(0, feature_string.length()-1);
		feature_wr.println(feature_final);
		feature_wr.flush();	
		
		
		
		
		
	}
	
	
}
