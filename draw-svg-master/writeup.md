#### Team member:

| Name        | SUNET ids |
| ----------- | --------- |
| Mingfei Guo | mfguo     |
| Boyu Zhang  | bzhang99  |

#### Tasks Breakdown:

Mingfei Guo: Task 1- 3

Boyu Zhang: Task 4 - 6

#### Tasks comments:

Mingfei Guo:

Boyu Zhang: 

For task 4, I first implemented the "sample_nearest" function to ensure the "rasterize_image" function I wrote was correct. Then, I follow the bilinear interpolation instruction on the slide to implement the "sample_bilinear" function. The tricky part is we need to subtract 0.5 from u and v to get the correct u00, u10, u01, and u11 coordinate.

For task 5, I followed the "Simple alpha compositing" instruction to implement the alpha_blending function. The tricky part is we need to pre-multiply the alpha value on both element color and canvas color.

For task 6, I used adobe illustrator to draw a tree. I had to modify the SVG file manually because adobe illustrator used another way to specify the polygon fill and stroke color.s
