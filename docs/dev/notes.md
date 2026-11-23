Starting to actually write down some of the notes so that I can keep track of them and they are accessible to LLMs and the like.

Material transparency:

- For EEVEE we have to choose between Blended and Dithered. Blended doesn't work with AO and transparency, no AO is displayed. Dithered works with AO but doesn't work with the depth-based outline method, as we get speckles of depth before is intermittently has transparency instead of block gradual.

"Transparency Overlap" needs to be disabled for `Flat` as otherwise introduces sorting errors in the transparency.
