# https://github.com/hexgrad/kokoro
# pip install -q kokoro>=0.7.16 soundfile
# pip install IPython
# pip index versions kokoro

from kokoro import KPipeline
from IPython.display import display, Audio
import soundfile as sf
import torch
import numpy as np

pipeline = KPipeline(lang_code="a")

from text_TTS import text

generator = pipeline(text, voice="af_heart", speed=1.0, split_pattern=r"\n+")

audio_segments = []

for i, (gs, ps, audio) in enumerate(generator):
    audio_segments.append(audio)

combined_audio = np.concatenate(audio_segments)

sf.write("combined_audio.wav", combined_audio, 24000)
