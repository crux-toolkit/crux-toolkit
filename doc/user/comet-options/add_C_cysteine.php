<?php include "head.php" ; ?>
<html>
<body>
<?php include "topmenu.php" ; ?>
<?php include "imgbar.php" ; ?>

<div id="page">
   <div id="content_full">
      <div class="post hr">

         <h2>Comet parameter: add_C_cysteine</h2>

         <ul>
         <li>Specify a static modification to the residue C.
         <li>The specified mass is added to the unmodified mass of C.
         <li>The default value is "0.0" if this parameter is missing.
         </ul>

         <p>Example:
         <br><tt>add_C_cysteine = 57.02</tt>

      </div>
   </div>
   <div style="clear: both;">&nbsp;</div>
</div>

<?php include "footer.php" ; ?>
</html>
