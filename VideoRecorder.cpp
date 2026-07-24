#include "VideoRecorder.h"
#include <iostream>
#include <cstdlib>

VideoRecorder::VideoRecorder(const char *filename, int width, int height, int fps)
	: fmtCtx(nullptr), codecCtx(nullptr), stream(nullptr),
	  frame(nullptr), pkt(nullptr), swsCtx(nullptr),
	  width(width), height(height), frameIndex(0)
{
	avformat_alloc_output_context2(&fmtCtx, nullptr, nullptr, filename);
	if (!fmtCtx) {
		std::cerr << "VideoRecorder: Could not allocate output context" << std::endl;
		std::exit(1);
	}

	const AVCodec *codec = avcodec_find_encoder(AV_CODEC_ID_H264);
	if (!codec) {
		std::cerr << "VideoRecorder: H.264 encoder not found" << std::endl;
		std::exit(1);
	}

	stream = avformat_new_stream(fmtCtx, nullptr);
	if (!stream) {
		std::cerr << "VideoRecorder: Could not create stream" << std::endl;
		std::exit(1);
	}

	codecCtx = avcodec_alloc_context3(codec);
	codecCtx->codec_id = AV_CODEC_ID_H264;
	codecCtx->width = width;
	codecCtx->height = height;
	codecCtx->time_base = {1, fps};
	codecCtx->framerate = {fps, 1};
	codecCtx->pix_fmt = AV_PIX_FMT_YUV420P;
	codecCtx->color_range = AVCOL_RANGE_JPEG;
	codecCtx->colorspace = AVCOL_SPC_BT709;
	codecCtx->color_primaries = AVCOL_PRI_BT709;
	codecCtx->color_trc = AVCOL_TRC_IEC61966_2_1;
	codecCtx->gop_size = 12;
	codecCtx->max_b_frames = 2;

	if (fmtCtx->oformat->flags & AVFMT_GLOBALHEADER)
		codecCtx->flags |= AV_CODEC_FLAG_GLOBAL_HEADER;

	av_opt_set(codecCtx->priv_data, "preset", "fast", 0);
	av_opt_set(codecCtx->priv_data, "crf", "23", 0);
	av_opt_set(codecCtx->priv_data, "x264-params",
		"colorprim=bt709:transfer=iec61966-2-1:colormatrix=bt709:fullrange=on",
		0);

	if (avcodec_open2(codecCtx, codec, nullptr) < 0) {
		std::cerr << "VideoRecorder: Could not open codec" << std::endl;
		std::exit(1);
	}

	// [REMOVED DUPLICATE HEADER CALL THAT CAUSES THE CODEC NONE ERROR FROM HERE]

	frame = av_frame_alloc();
	frame->format = codecCtx->pix_fmt;
	frame->width = width;
	frame->height = height;
	frame->color_range = AVCOL_RANGE_JPEG;
	frame->colorspace = AVCOL_SPC_BT709;
	frame->color_primaries = AVCOL_PRI_BT709;
	frame->color_trc = AVCOL_TRC_IEC61966_2_1;
	av_frame_get_buffer(frame, 0);

	pkt = av_packet_alloc();

	swsCtx = sws_getContext(
		width, height, AV_PIX_FMT_RGB24,
		width, height, AV_PIX_FMT_YUV420P,
		SWS_BILINEAR, nullptr, nullptr, nullptr
	);

	int *inv_table, *table;
	int srcRange, dstRange, brightness, contrast, saturation;

	if (sws_getColorspaceDetails(swsCtx, &inv_table, &srcRange, &table, &dstRange,
								 &brightness, &contrast, &saturation) >= 0)
	{
		srcRange = 1; // PC/Full range for source RGB
		dstRange = 1; // PC/Full range for output YUV

		sws_setColorspaceDetails(
			swsCtx,
			sws_getCoefficients(SWS_CS_ITU709), srcRange,
			sws_getCoefficients(SWS_CS_ITU709), dstRange,
			brightness, contrast, saturation
		);
	}

	// Copy basic encoder specs over to the stream structure
	avcodec_parameters_from_context(stream->codecpar, codecCtx);
	
	// Enforce the codec parameters explicitly so the MP4 muxer maps the stream
	stream->codecpar->codec_type      = AVMEDIA_TYPE_VIDEO;
	stream->codecpar->codec_id        = codecCtx->codec_id;

	// Push color space specifications to the container parameters
	stream->codecpar->color_range     = codecCtx->color_range;
	stream->codecpar->color_space     = codecCtx->colorspace;
	stream->codecpar->color_primaries = codecCtx->color_primaries;
	stream->codecpar->color_trc       = codecCtx->color_trc;

	stream->time_base = codecCtx->time_base;

	// Open file and write header ONLY ONCE after parameters are fully bound
	if (!(fmtCtx->oformat->flags & AVFMT_NOFILE)) {
		if (avio_open(&fmtCtx->pb, filename, AVIO_FLAG_WRITE) < 0) {
			std::cerr << "VideoRecorder: Could not open output file" << std::endl;
			std::exit(1);
		}
	}

	if (avformat_write_header(fmtCtx, nullptr) < 0) {
		std::cerr << "VideoRecorder: Could not write header" << std::endl;
		std::exit(1);
	}
}


VideoRecorder::~VideoRecorder()
{
	avcodec_send_frame(codecCtx, nullptr);
	WritePackets();

	av_write_trailer(fmtCtx);

	if (!(fmtCtx->oformat->flags & AVFMT_NOFILE))
		avio_closep(&fmtCtx->pb);

	sws_freeContext(swsCtx);
	av_frame_free(&frame);
	av_packet_free(&pkt);
	avcodec_free_context(&codecCtx);
	avformat_free_context(fmtCtx);
}

void VideoRecorder::WritePackets()
{
	while (avcodec_receive_packet(codecCtx, pkt) == 0) {
		av_packet_rescale_ts(pkt, codecCtx->time_base, stream->time_base);
		pkt->stream_index = stream->index;
		av_interleaved_write_frame(fmtCtx, pkt);
		av_packet_unref(pkt);
	}
}

void VideoRecorder::WriteFrame(const uint8_t *rgbData)
{
	av_frame_make_writable(frame);

	const uint8_t *srcSlice[1] = { rgbData };
	int srcStride[1] = { 3 * width };

	sws_scale(swsCtx, srcSlice, srcStride, 0, height,
			  frame->data, frame->linesize);

	frame->pts = frameIndex++;

	avcodec_send_frame(codecCtx, frame);
	WritePackets();
}
