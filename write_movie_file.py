import os
import moviepy.video.io.ImageSequenceClip
image_folder='/home/d9smith/projects/age/volume/plots/test_images'
fps=1

image_files = [os.path.join(image_folder,img)
               for img in os.listdir(image_folder)
               if img.endswith(".png")]
clip = moviepy.video.io.ImageSequenceClip.ImageSequenceClip(image_files, fps=fps)
clip.write_videofile('/home/d9smith/projects/age/volume/my_test_video.mp4')