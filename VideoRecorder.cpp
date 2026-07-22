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
	codecCtx->gop_size = 12;
	codecCtx->max_b_frames = 2;

	if (fmtCtx->oformat->flags & AVFMT_GLOBALHEADER)
		codecCtx->flags |= AV_CODEC_FLAG_GLOBAL_HEADER;

	av_opt_set(codecCtx->priv_data, "preset", "fast", 0);
	av_opt_set(codecCtx->priv_data, "crf", "23", 0);

	if (avcodec_open2(codecCtx, codec, nullptr) < 0) {
		std::cerr << "VideoRecorder: Could not open codec" << std::endl;
		std::exit(1);
	}

	avcodec_parameters_from_context(stream->codecpar, codecCtx);
	stream->time_base = codecCtx->time_base;

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

	frame = av_frame_alloc();
	frame->format = codecCtx->pix_fmt;
	frame->width = width;
	frame->height = height;
	av_frame_get_buffer(frame, 0);

	pkt = av_packet_alloc();

	swsCtx = sws_getContext(
		width, height, AV_PIX_FMT_RGB24,
		width, height, AV_PIX_FMT_YUV420P,
		SWS_BILINEAR, nullptr, nullptr, nullptr
	);
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
