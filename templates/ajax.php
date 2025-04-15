<?php
    $con=mysqli_connect('localhost','root','','demo')
    if(isset($_REQUEST['hdn_insert']))
    {
        $data=array();
        $input_1=$_REQUEST['input1_pl'];
        $input_2=$_REQUEST['input2_pl'];

        $sel="insert into info1 values(NULL,'$input_1', '$input_2')";
        mysqli_query($con,$sel);

        if(mysqli_affected_rows($con))
        {
            $data['status']=1;
        }
        else
        {
            $data['status']=0;
        }
        echo json_encode($data)
    }
?>