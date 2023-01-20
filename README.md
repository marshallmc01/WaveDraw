# WaveDraw
Personal Project I'm working on :)

# description

The idea is that you can draw a wave/line/squiggle to form a base. With that base, you can manipulate it into your ideal waveform via some functions I've added. Once the wave is finished, you can save it as a vector, play it as audio, superimpose it with other waves, etc.

# Highlights

The main function is a mess so good luck navigating it! I recommend you start by looking at the regressions.cpp file, it contains some super cool algorithms for bezier regressions. After that, the split/merge fragments logic is pretty cool. It's not perfect but it works and I think some parts are clever.

# status and future plan

Currently I've finished basic drawing features (though there's so much more I could add).
Moving forwards, I'm going to create a new window for superimposing waves, as well as messing with a wave's frequency and amplitude.
After that I'll get to the audio (which I may mess around with earlier rather than later for fun).

# Repo navigation

I'm new to building actual projects so right now I'm just leaving an exe file and the source code. In the future I'll make it more professional but for now i'm just having fun.

# Credits

I used a fair amount of resources to develop the regressions, so I've listed what I found helpful below.
&nbsp;  
&nbsp;  

Bez regression in C: https://github.com/erich666/GraphicsGems/blob/master/gems/FitCurves.c  
Constrained least squares (polyfit with fixed endpoints): http://www.seas.ucla.edu/~vandenbe/133A/lectures/cls.pdf  
Matrix solution for cubic bezier regression: https://www.jimherold.com/computer-science/best-fit-bezier-curve   
Visualization of newton raphson reparamaterization: https://stackoverflow.com/questions/42315048/bezier-curve-fitting-with-known-end-points/42325165#42325165  
Polynomial regression maths: https://muthu.co/maths-behind-polynomial-regression/  





