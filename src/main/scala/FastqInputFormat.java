import java.io.IOException;
import java.util.Vector;

import org.apache.hadoop.classification.InterfaceAudience;
import org.apache.hadoop.classification.InterfaceStability;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FSDataInputStream;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.fs.Seekable;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.io.compress.CodecPool;
import org.apache.hadoop.io.compress.CompressionCodec;
import org.apache.hadoop.io.compress.CompressionCodecFactory;
import org.apache.hadoop.io.compress.Decompressor;
import org.apache.hadoop.io.compress.SplitCompressionInputStream;
import org.apache.hadoop.io.compress.SplittableCompressionCodec;
import org.apache.hadoop.mapreduce.InputSplit;
import org.apache.hadoop.mapreduce.RecordReader;
import org.apache.hadoop.mapreduce.TaskAttemptContext;
import org.apache.hadoop.mapreduce.lib.input.FileSplit;
import org.apache.hadoop.mapreduce.lib.input.TextInputFormat;
import org.apache.hadoop.util.LineReader;

/*
 * Description: 
 * 		A new API hadoop InputFormat for Fastq file loading.
 * Usage:
import org.apache.hadoop.conf.Configuration
import org.apache.hadoop.io.{LongWritable, Text}
import org.apache.hadoop.mapreduce.lib.input.TextInputFormat
import org.apache.spark.{SparkContext, SparkConf} // for shell
val conf = new SparkConf().setAppName("decos").registerKryoClasses(
	Array[Class[_]](Class.forName("org.apache.hadoop.io.LongWritable"),
					Class.forName("org.apache.hadoop.io.Text")
					))
val dataset = sc.newAPIHadoopFile("/home/x/xieluyu/reads/Ecoli1.fq", classOf[FastqInputFormat], classOf[LongWritable], classOf[Text])
val data = dataset.map(x=>(x._1.get,x._2.toString))
data.count()
data.collect()


 * By: 
 * 		Luyu
 */

public class FastqInputFormat extends TextInputFormat {
	public RecordReader<LongWritable, Text> createRecordReader(InputSplit split, TaskAttemptContext context) {
		return new FastqRecordReader();
	}
}

@InterfaceAudience.LimitedPrivate({ "MapReduce", "Pig" })
@InterfaceStability.Evolving
class FastqRecordReader extends RecordReader<LongWritable, Text> {
	public static final String MAX_LINE_LENGTH = "mapreduce.input.linerecordreader.line.maxlength";

	private long start;
	private long pos;
	private long end;
	private LineReader in;
	private FSDataInputStream fileIn;
	private Seekable filePosition;
	private int maxLineLength;
	private LongWritable key;
	private Text value;
	private boolean isCompressedInput;
	private Decompressor decompressor;
	private byte[] recordDelimiterBytes;
	private Vector<String> preread = new Vector<String>();
	private Vector<Integer> prereadlen = new Vector<Integer>();

	public FastqRecordReader() {
	}

	public FastqRecordReader(byte[] recordDelimiter) {
		this.recordDelimiterBytes = recordDelimiter;
	}

	public void initialize(InputSplit genericSplit, TaskAttemptContext context) throws IOException {
		FileSplit split = (FileSplit) genericSplit;
		Configuration job = context.getConfiguration();
		this.maxLineLength = job.getInt(MAX_LINE_LENGTH, Integer.MAX_VALUE);
		start = split.getStart();
		end = start + split.getLength();
		final Path file = split.getPath();

		// open the file and seek to the start of the split
		final FileSystem fs = file.getFileSystem(job);
		fileIn = fs.open(file);

		CompressionCodec codec = new CompressionCodecFactory(job).getCodec(file);
		if (null != codec) {
			isCompressedInput = true;
			decompressor = CodecPool.getDecompressor(codec);
			if (codec instanceof SplittableCompressionCodec) {
				final SplitCompressionInputStream cIn = ((SplittableCompressionCodec) codec).createInputStream(
						fileIn, decompressor, start, end,
						SplittableCompressionCodec.READ_MODE.BYBLOCK);
				if (null == this.recordDelimiterBytes) {
					in = new LineReader(cIn, job);
				} else {
					in = new LineReader(cIn, job, this.recordDelimiterBytes);
				}

				start = cIn.getAdjustedStart();
				end = cIn.getAdjustedEnd();
				filePosition = cIn;
			} else {
				if (null == this.recordDelimiterBytes) {
					in = new LineReader(codec.createInputStream(fileIn, decompressor),
							job);
				} else {
					in = new LineReader(codec.createInputStream(fileIn,
							decompressor), job, this.recordDelimiterBytes);
				}
				filePosition = fileIn;
			}
		} else {
			fileIn.seek(start);
			if (null == this.recordDelimiterBytes) {
				in = new LineReader(fileIn, job);
			} else {
				in = new LineReader(fileIn, job, this.recordDelimiterBytes);
			}
			filePosition = fileIn;
		}
		// If this is not the first split, we always throw away first record
		// because we always (except the last split) read one extra line in
		// next() method.
		if (start != 0) {
			start += readChunkedFastq(maxBytesToConsume(start));
		}
		this.pos = start;
	}

	private int readChunkedFastq(int maxBytesToConsume) throws IOException {
		Text line = new Text();
		int len = 0;
		while (true) {
			len = in.readLine(line, this.maxLineLength, maxBytesToConsume);
			maxBytesToConsume -= len;
			if (len == 0) return 0;
			preread.add(line.toString());
			prereadlen.add(len);
			if (preread.size() >3 && line.getLength()>0 && line.charAt(0) == '+' &&
                    preread.get(preread.size()-3).length()>0 && preread.get(preread.size()-3).charAt(0)=='@') break;
		}
//		System.out.println("[" + start + "," + end + "] buffer size :" + preread.size());
//		for (Text i : preread)
//			System.out.println("[" + start + "," + end + "]:" + i.toString());
		for (len = 0; preread.size() > 3;) {
			len += prereadlen.remove(0);
			preread.remove(0);
		}
		return len;
	}

	private int readFastq(Text str, int maxLineLength, int maxBytesToConsume) throws IOException {
		Text line = new Text();
		StringBuilder fastq = new StringBuilder();
		int len = 0, fastqLen = 0;
		for (int i = 0; i < 4; i++) {
			if (preread.isEmpty()) {
                len = in.readLine(line, this.maxLineLength, maxBytesToConsume);
                fastq.append(line.toString()).append('\n');
            } else {
				fastq.append(preread.remove(0)).append('\n');
				len = prereadlen.remove(0);
			}
			maxBytesToConsume -= len;
			fastqLen += len;
		}
		str.set(fastq.toString());
		return fastqLen;
	}

	private int maxBytesToConsume(long pos) {
		return isCompressedInput
				? Integer.MAX_VALUE
				: (int) Math.min(Integer.MAX_VALUE, end - pos);
	}

	private long getFilePosition() throws IOException {
		long retVal;
		if (isCompressedInput && null != filePosition) {
			retVal = filePosition.getPos();
		} else {
			retVal = pos;
		}
		return retVal;
	}

	public boolean nextKeyValue() throws IOException {
		if (key == null) {
			key = new LongWritable();
		}
		key.set(pos);
		if (value == null) {
			value = new Text();
		}
		int newSize = 0;
		// We always read one extra line, which lies outside the upper
		// split limit i.e. (end - 1)
		while (getFilePosition() <= end) {
			newSize = readFastq(value, maxLineLength, Math.max(maxBytesToConsume(pos), maxLineLength));
			pos += newSize;
			if (newSize < maxLineLength) {
				break;
			}
		}
		if (newSize == 0) {
			key = null;
			value = null;
			return false;
		} else {
			return true;
		}
	}

	@Override
	public LongWritable getCurrentKey() {
		return key;
	}

	@Override
	public Text getCurrentValue() {
		return value;
	}

	/**
	 * Get the progress within the split
	 */
	public float getProgress() throws IOException {
		if (start == end) {
			return 0.0f;
		} else {
			return Math.min(1.0f, (getFilePosition() - start) / (float) (end - start));
		}
	}

	public synchronized void close() throws IOException {
		try {
			if (in != null) {
				in.close();
			}
		} finally {
			if (decompressor != null) {
				CodecPool.returnDecompressor(decompressor);
			}
		}
	}
}
