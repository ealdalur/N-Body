#pragma once

#include <cstdint>

extern "C" {
#include <libavcodec/avcodec.h>
#include <libavformat/avformat.h>
#include <libavutil/imgutils.h>
#include <libavutil/opt.h>
#include <libswscale/swscale.h>
}

class VideoRecorder
{
	AVFormatContext *fmtCtx;
	AVCodecContext *codecCtx;
	AVStream *stream;
	AVFrame *frame;
	AVPacket *pkt;
	SwsContext *swsCtx;

	int width, height;
	int64_t frameIndex;

	void WritePackets();

public:
	VideoRecorder(const char *filename, int width, int height, int fps);
	~VideoRecorder();

	void WriteFrame(const uint8_t *rgbData);
};
